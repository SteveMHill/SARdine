use std::collections::HashMap;

use crate::io::annotation::AnnotationRoot;
use crate::types::{BurstTiming, SarResult, StateVector, SubSwath};

/// Doppler centroid polynomial model derived from annotation metadata
#[derive(Debug, Clone)]
pub struct DopplerCentroidModel {
    /// Reference time (seconds since start) for the polynomial origin
    pub t0: f64,
    /// Polynomial coefficients ordered from constant term upwards
    pub coeffs: Vec<f64>,
}

impl DopplerCentroidModel {
    /// Evaluate the Doppler centroid polynomial at a given azimuth time (seconds since product start)
    pub fn evaluate(&self, azimuth_time: f64) -> f64 {
        let mut acc = 0.0;
        let mut pow = 1.0;
        let x = azimuth_time - self.t0;
        for coeff in &self.coeffs {
            acc += coeff * pow;
            pow *= x;
        }
        acc
    }
}

/// Pre-computed burst timing segment for fast azimuth mapping
#[derive(Debug, Clone)]
pub struct BurstSegment {
    pub subswath_id: String,
    pub burst_index: usize,
    pub start_time_rel: f64,
    pub end_time_rel: f64,
    pub start_line: f64,
    pub end_line: f64,
    pub line_time_interval: f64,
}

/// Range-Doppler terrain correction parameters
#[derive(Debug, Clone)]
pub struct RangeDopplerParams {
    /// NATIVE range pixel spacing in meters (sensor physics, NOT multilooked grid spacing)
    /// For Sentinel-1 IW: ~2.33m (native) not ~9.3m (4x multilooked)
    pub range_pixel_spacing: f64,
    /// NATIVE azimuth time interval in seconds (sensor physics, NOT grid spacing in meters)
    /// For Sentinel-1 TOPS: subswath-specific PRF interval ~0.0006s, not grid spacing ~14m
    pub azimuth_pixel_spacing: f64,
    /// Slant range time to first pixel (seconds)
    pub slant_range_time: f64,
    /// Pulse repetition frequency (Hz)
    pub prf: f64,
    /// Azimuth time interval (seconds per line) - CRITICAL for TOPS data
    /// This is the actual line time from annotation, NOT necessarily 1/PRF
    pub azimuth_time_interval: f64,
    /// Radar wavelength (meters)
    pub wavelength: f64,
    /// Speed of light (m/s)
    pub speed_of_light: f64,
    /// CRITICAL: Orbit reference epoch in seconds since Unix epoch
    /// All times must be computed relative to this to avoid solver failures
    pub orbit_ref_epoch_utc: f64,
    /// Product start time in seconds since orbit_ref_epoch (NOT Unix epoch)
    pub product_start_rel_s: f64,
    /// Absolute product start time in seconds (Unix epoch) - for legacy compatibility
    #[deprecated(
        since = "0.2.2",
        note = "Use product_start_rel_s with orbit_ref_epoch_utc"
    )]
    pub product_start_time_abs: f64,
    /// Absolute product stop time in seconds (Unix epoch) - for legacy compatibility
    #[deprecated(since = "0.2.2", note = "Use product_start_rel_s + product_duration")]
    pub product_stop_time_abs: f64,
    /// Total acquisition duration in seconds (stop - start)
    pub product_duration: f64,
    /// Total azimuth lines in the (merged) product grid (metadata authoritative)
    pub total_azimuth_lines: Option<usize>,
    /// Optional Doppler centroid polynomial model
    pub doppler_centroid: Option<DopplerCentroidModel>,
    /// First valid line (azimuth) - TOPS burst ramp-up exclusion
    pub first_valid_line: Option<usize>,
    /// Last valid line (azimuth) - TOPS burst ramp-down exclusion
    pub last_valid_line: Option<usize>,
    /// First valid sample (range) - swath boundary
    pub first_valid_sample: Option<usize>,
    /// Last valid sample (range) - swath boundary
    pub last_valid_sample: Option<usize>,
    /// Range multilook factor (native pixels per output pixel)
    pub range_multilook_factor: f64,
    /// Azimuth multilook factor (native lines per output line)
    pub azimuth_multilook_factor: f64,
    /// Pre-computed safe range multilook factor (always >= 1.0)
    /// This avoids repeated `.max(1.0)` calculations in hot loops
    pub range_multilook_safe: f64,
    /// Pre-computed safe azimuth multilook factor (always >= 1.0)
    /// This avoids repeated `.max(1.0)` calculations in hot loops
    pub azimuth_multilook_safe: f64,
    /// Subswath geometry information for per-subswath Range-Doppler calculation
    pub subswaths: HashMap<String, SubSwath>,
    /// Detailed burst timing table used for azimuth mapping
    pub burst_timings: Vec<BurstTiming>,
    /// Piecewise-linear segments derived from burst timings (multilooked domain)
    pub burst_segments: Vec<BurstSegment>,
    /// Reference (ellipsoid) incidence angle in degrees for RTC normalization
    /// Used by Small 2011 area-projection formula: γ⁰ = σ⁰ × cos(θ_ref) / (cos(θ_local) × M)
    /// If None, defaults to 35° (typical Sentinel-1 IW mid-swath)
    pub reference_incidence_angle_deg: Option<f64>,
    /// Near-range (first sample) ellipsoid incidence angle in degrees
    /// For Sentinel-1 IW: ~29.1° (IW1) to ~39.3° (IW3 near edge)
    pub incidence_angle_near_deg: Option<f64>,
    /// Far-range (last sample) ellipsoid incidence angle in degrees  
    /// For Sentinel-1 IW: ~35.5° (IW1) to ~46.0° (IW3 far edge)
    pub incidence_angle_far_deg: Option<f64>,
    /// Total range samples for incidence angle interpolation
    pub total_range_samples: Option<usize>,
}

impl RangeDopplerParams {
    /// Compute and set safe multilook factors (always >= 1.0)
    /// This should be called after constructing RangeDopplerParams to pre-compute
    /// values that are used repeatedly in hot loops
    #[inline]
    pub fn compute_safe_multilook_factors(&mut self) {
        self.range_multilook_safe = self.range_multilook_factor.max(1.0);
        self.azimuth_multilook_safe = self.azimuth_multilook_factor.max(1.0);
    }

    /// Compute boundary tolerance based on actual subswath spacing.
    /// Returns tolerance in pixels - defaults to 20 if no subswaths available.
    ///
    /// The tolerance is set to ~5% of the minimum subswath spacing, clamped to [10, 50] pixels.
    /// This adapts to actual product geometry instead of using hardcoded values.
    pub fn compute_boundary_tolerance(&self) -> usize {
        const DEFAULT_TOLERANCE: usize = 20;
        const MIN_TOLERANCE: usize = 10;
        const MAX_TOLERANCE: usize = 50;

        if self.subswaths.len() < 2 {
            return DEFAULT_TOLERANCE;
        }

        // Collect subswath boundaries sorted by range
        let mut boundaries: Vec<(usize, usize)> = self
            .subswaths
            .values()
            .map(|s| (s.first_sample_global, s.last_sample_global))
            .collect();
        boundaries.sort_by_key(|(start, _)| *start);

        // Find minimum gap between adjacent subswaths
        let mut min_gap = usize::MAX;
        for window in boundaries.windows(2) {
            let gap = window[1].0.saturating_sub(window[0].1);
            if gap > 0 && gap < min_gap {
                min_gap = gap;
            }
        }

        if min_gap == usize::MAX {
            return DEFAULT_TOLERANCE;
        }

        // Use ~5% of minimum gap as tolerance, clamped to reasonable bounds
        let computed = (min_gap / 20).max(MIN_TOLERANCE).min(MAX_TOLERANCE);
        log::debug!(
            "Computed boundary tolerance: {} pixels (from min subswath gap of {} pixels)",
            computed,
            min_gap
        );
        computed
    }

    /// Absolute product start time in UTC seconds
    #[inline]
    pub fn product_start_absolute(&self) -> f64 {
        self.orbit_ref_epoch_utc + self.product_start_rel_s
    }

    /// Absolute product stop time in UTC seconds
    #[inline]
    pub fn product_stop_absolute(&self) -> f64 {
        self.product_start_absolute() + self.product_duration
    }

    /// Compute per-pixel ellipsoid incidence angle (cosine) from range sample position.
    ///
    /// This interpolates linearly between incidence_angle_near and incidence_angle_far
    /// based on the range sample index, providing range-dependent reference incidence
    /// for more accurate RTC normalization.
    ///
    /// # Arguments
    /// * `range_sample` - Range sample index (0-based, in native pixel coordinates)
    ///
    /// # Returns
    /// * Cosine of ellipsoid incidence angle at this range position
    ///
    /// # Scientific Rationale
    /// Using a constant mid-swath reference angle introduces ~1-2 dB range-dependent
    /// bias across the IW swath (29°-46° variation). Linear interpolation reduces
    /// this to sub-0.5 dB for typical Sentinel-1 acquisitions.
    #[inline]
    pub fn cos_ellipsoid_incidence_at_range(&self, range_sample: usize) -> f64 {
        // Check if we have the near/far angles and total samples
        let (near_deg, far_deg, total_samples) = match (
            self.incidence_angle_near_deg,
            self.incidence_angle_far_deg,
            self.total_range_samples,
        ) {
            (Some(near), Some(far), Some(total)) if total > 0 => (near, far, total),
            _ => {
                // Fallback to constant reference angle
                let ref_deg = self.reference_incidence_angle_deg.unwrap_or(35.0);
                return ref_deg.to_radians().cos();
            }
        };

        // Linear interpolation between near and far incidence angles
        let t = (range_sample as f64) / ((total_samples - 1).max(1) as f64);
        let t = t.clamp(0.0, 1.0);

        let incidence_deg = near_deg + t * (far_deg - near_deg);
        incidence_deg.to_radians().cos()
    }

    /// Build burst segments from raw timing metadata for fast azimuth lookup.
    ///
    /// # Critical: azimuth_time_rel_orbit Validation
    /// Each burst MUST have a valid (non-zero) `azimuth_time_rel_orbit` for accurate geocoding.
    /// If this value is 0.0, it typically indicates the JSON deserialization fell back to default,
    /// meaning the source `BurstRecord` had `null` for this field. This WILL cause wrong geocoding!
    /// Build burst segments from raw timing metadata for fast azimuth lookup.
    ///
    /// # Parameters
    /// - `native_azimuth_dt`: Within-burst line spacing in seconds (= 1/PRF).
    ///   NOT the product-average `azimuth_time_interval` which includes inter-burst gaps!
    /// - `total_debursted_lines`: Total lines in post-deburst image. Used to compute
    ///   correct cumulative line offsets. If None, falls back to annotation line counts.
    ///
    /// For merged TOPSAR, burst segments are built PER SUBSWATH because IW1/IW2/IW3
    /// are acquired at different (staggered) times. Each subswath's bursts independently
    /// map time→line, but all map to the same output image lines.
    pub fn build_burst_segments(
        bursts: &[BurstTiming],
        native_azimuth_dt: f64,
        azimuth_multilook_factor: f64,
        total_debursted_lines: Option<usize>,
    ) -> Vec<BurstSegment> {
        if bursts.is_empty() || !native_azimuth_dt.is_finite() || native_azimuth_dt <= 0.0 {
            return Vec::new();
        }

        // CRITICAL VALIDATION: Check for zero azimuth_time_rel_orbit values
        // Note: First burst at exactly 0.0 is legitimate (orbit-relative timing starts at epoch)
        // We only fail if ALL bursts are zero (indicates deserialization bug)
        let zero_count = bursts
            .iter()
            .filter(|b| b.azimuth_time_rel_orbit.abs() < 1e-9)
            .count();

        if zero_count == bursts.len() && bursts.len() > 0 {
            log::error!(
                "❌ CRITICAL: ALL {} bursts have azimuth_time_rel_orbit ≈ 0.0! (deserialization bug)",
                bursts.len()
            );
            for (idx, burst) in bursts.iter().take(5).enumerate() {
                log::error!(
                    "   Burst #{}: subswath={}, burst_index={}, azimuth_time_rel_orbit={:.9}",
                    idx + 1, burst.subswath_id, burst.burst_index, burst.azimuth_time_rel_orbit
                );
            }
            log::error!("   ❌ FAIL-FAST: Returning empty burst segments.");
            return Vec::new();
        } else if zero_count > 0 {
            log::debug!(
                "📐 Found {} of {} bursts with azimuth_time_rel_orbit ≈ 0.0 (first burst or epoch-relative)",
                zero_count,
                bursts.len()
            );
        }

        let ml_factor = if azimuth_multilook_factor.is_finite() && azimuth_multilook_factor > 0.0 {
            azimuth_multilook_factor
        } else {
            1.0
        };

        // Deduplicate by (subswath_id, burst_index, first_line_global) to remove polarization duplicates
        let mut seen = std::collections::HashSet::new();
        let deduped_bursts: Vec<&BurstTiming> = bursts
            .iter()
            .filter(|b| {
                let key = (b.subswath_id.clone(), b.burst_index, b.first_line_global);
                seen.insert(key)
            })
            .collect();

        if deduped_bursts.len() < bursts.len() {
            log::warn!(
                "📐 Deduplicated burst records: {} → {} (removed {} duplicates)",
                bursts.len(),
                deduped_bursts.len(),
                bursts.len() - deduped_bursts.len()
            );
        }

        // Group bursts by subswath to build per-subswath segments
        let mut by_subswath: std::collections::BTreeMap<String, Vec<&BurstTiming>> =
            std::collections::BTreeMap::new();
        for b in &deduped_bursts {
            by_subswath
                .entry(b.subswath_id.clone())
                .or_default()
                .push(b);
        }

        let num_subswaths = by_subswath.len();

        // Log subswath grouping for debugging
        let sw_counts: Vec<String> = by_subswath
            .iter()
            .map(|(k, v)| format!("{}:{}", k, v.len()))
            .collect();
        log::warn!(
            "📐 Burst segments: {} unique records, {} subswaths [{}], native_azimuth_dt={:.9}s",
            deduped_bursts.len(),
            num_subswaths,
            sw_counts.join(", "),
            native_azimuth_dt,
        );
        let mut all_segments: Vec<BurstSegment> = Vec::new();

        for (sw_id, sw_bursts) in &mut by_subswath {
            // Sort by azimuth time within this subswath
            sw_bursts.sort_by(|a, b| {
                a.azimuth_time_rel_orbit
                    .partial_cmp(&b.azimuth_time_rel_orbit)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

            let num_bursts = sw_bursts.len();

            // Each subswath independently maps time→line, but all produce the same number
            // of output lines. valid_lines_per_burst is per-subswath's number of bursts.
            let valid_lines_per_burst: f64 = if let Some(total) = total_debursted_lines {
                total as f64 / num_bursts as f64
            } else {
                sw_bursts
                    .iter()
                    .map(|b| {
                        (b.last_line_global.max(b.first_line_global) as f64
                            - b.first_line_global as f64)
                            .max(0.0)
                    })
                    .sum::<f64>()
                    / num_bursts as f64
            };

            // Compute the average annotation burst span to determine edge margin
            let avg_annotation_span: f64 = sw_bursts
                .iter()
                .map(|b| {
                    (b.last_line_global.max(b.first_line_global) as f64
                        - b.first_line_global as f64)
                        .max(0.0)
                })
                .sum::<f64>()
                / num_bursts as f64;

            // Margin: lines trimmed from each edge during deburst
            let margin_lines = ((avg_annotation_span - valid_lines_per_burst) / 2.0).max(0.0);
            let margin_time = margin_lines * native_azimuth_dt;

            log::warn!(
                "📐 Subswath {}: {} bursts, {:.1} valid lines/burst (of {:.0} annotation span), margin={:.0} lines, total_debursted_lines={:?}",
                sw_id, num_bursts, valid_lines_per_burst, avg_annotation_span, margin_lines, total_debursted_lines
            );

            // Build segments with POST-DEBURST cumulative line numbers
            let mut cumulative_lines: f64 = 0.0;
            for burst in sw_bursts {
                let valid_start_time = burst.azimuth_time_rel_orbit + margin_time;
                let valid_duration = valid_lines_per_burst * native_azimuth_dt;

                let start_line = cumulative_lines / ml_factor;
                let end_line = (cumulative_lines + valid_lines_per_burst) / ml_factor;
                let line_time_interval = native_azimuth_dt * ml_factor;

                cumulative_lines += valid_lines_per_burst;

                all_segments.push(BurstSegment {
                    subswath_id: burst.subswath_id.clone(),
                    burst_index: burst.burst_index,
                    start_time_rel: valid_start_time,
                    end_time_rel: valid_start_time + valid_duration,
                    start_line,
                    end_line,
                    line_time_interval,
                });
            }
        }

        all_segments.sort_by(|a, b| {
            a.start_time_rel
                .partial_cmp(&b.start_time_rel)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        all_segments
    }

    /// Derive an effective azimuth time interval that matches the actually observed image height.
    ///
    /// For merged TOPS data the metadata total_azimuth_lines can disagree with the debursted
    /// array height by a few percent due to burst-gap removal or trimming. When that happens we
    /// prefer a duration/line interval based on the observed multilooked height to keep
    /// azimuth-coordinate mapping inside bounds.
    pub fn effective_azimuth_time_interval(
        &self,
        observed_multilooked_lines: Option<usize>,
    ) -> f64 {
        // Use pre-computed safe value
        let az_factor = self.azimuth_multilook_safe;
        let observed_native = observed_multilooked_lines
            .map(|h| (h as f64) * az_factor)
            .filter(|v| *v > 0.0);
        let metadata_native = self.total_azimuth_lines.map(|v| v as f64);

        let candidate_native = match (metadata_native, observed_native) {
            (Some(meta), Some(obs)) => {
                let rel = ((meta - obs) / meta).abs();
                if rel > 0.02 {
                    log::warn!(
                        "⚠️ Azimuth line mismatch: metadata={} vs observed_ml={} (native≈{:.1}); using observed to derive effective line interval",
                        meta,
                        observed_multilooked_lines.unwrap_or(0),
                        obs
                    );
                    Some(obs)
                } else {
                    Some(meta)
                }
            }
            (Some(meta), None) => Some(meta),
            (None, Some(obs)) => Some(obs),
            (None, None) => None,
        };

        if let Some(native_lines) = candidate_native {
            if self.product_duration.is_finite()
                && self.product_duration > 0.0
                && native_lines > 0.0
            {
                return self.product_duration / native_lines;
            }
        }

        self.azimuth_time_interval
    }

    /// Create parameters from real annotation data with orbit vectors
    /// CRITICAL: Requires orbit vectors to establish orbit_ref_epoch
    /// This prevents accidental use of hardcoded parameters AND ensures proper time base
    pub fn from_annotation(
        annotation: &AnnotationRoot,
        orbit_vectors: &[StateVector],
    ) -> SarResult<Self> {
        annotation.extract_range_doppler_params(orbit_vectors)
    }

    /// Find the subswath that contains a given native range sample
    ///
    /// For merged IW data (IW1+IW2+IW3), returns the SubSwath whose sample range
    /// contains the given native_range_sample. This is critical for using the
    /// correct slant_range_time per subswath.
    ///
    /// Returns None if no subswath contains the sample (gap or outside bounds).
    #[inline]
    pub fn find_subswath_for_sample(&self, native_range_sample: f64) -> Option<&SubSwath> {
        let sample = native_range_sample as usize;

        // Subswaths are sorted by first_sample_global in IW mode:
        // IW1: 0..~22k, IW2: ~22k..~45k, IW3: ~45k..~67k
        // CRITICAL FIX: last_sample_global is EXCLUSIVE (one past the end), so use < not <=
        for subswath in self.subswaths.values() {
            if sample >= subswath.first_sample_global && sample < subswath.last_sample_global {
                return Some(subswath);
            }
        }

        // IMPROVED EDGE CASE HANDLING: Use computed boundary tolerance based on actual subswath spacing
        // This handles:
        // - Overlap regions between subswaths (IW1/IW2, IW2/IW3)
        // - Edge pixels that fall just outside bounds due to geoid/DEM errors
        // - Multilooking effects that may shift boundaries slightly
        let boundary_tolerance = self.compute_boundary_tolerance();
        let sample_f64 = native_range_sample;

        // Check if sample is within tolerance of any subswath boundary
        for subswath in self.subswaths.values() {
            let dist_to_start = if sample_f64 < subswath.first_sample_global as f64 {
                subswath.first_sample_global as f64 - sample_f64
            } else {
                f64::INFINITY
            };
            let dist_to_end = if sample_f64 >= subswath.last_sample_global as f64 {
                sample_f64 - (subswath.last_sample_global as f64)
            } else {
                f64::INFINITY
            };

            let min_dist = dist_to_start.min(dist_to_end);
            if min_dist <= boundary_tolerance as f64 {
                // Sample is within tolerance of this subswath's boundary
                static BOUNDARY_FALLBACK_COUNT: std::sync::atomic::AtomicUsize =
                    std::sync::atomic::AtomicUsize::new(0);
                let count =
                    BOUNDARY_FALLBACK_COUNT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if count < 20 {
                    log::debug!(
                        "⚠️  Sample {} near {} boundary (dist={:.1} pixels), using subswath {}",
                        sample,
                        if dist_to_start < dist_to_end {
                            "start"
                        } else {
                            "end"
                        },
                        min_dist,
                        subswath.id
                    );
                }
                return Some(subswath);
            }
        }

        // DIAGNOSTIC: Log when sample doesn't match any subswath (potential gap or boundary issue)
        static MISSING_SUBSWATH_LOG_COUNT: std::sync::atomic::AtomicUsize =
            std::sync::atomic::AtomicUsize::new(0);
        let log_count =
            MISSING_SUBSWATH_LOG_COUNT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if log_count < 20 {
            let mut subswath_ranges = String::new();
            for (name, sw) in &self.subswaths {
                subswath_ranges.push_str(&format!(
                    "{}:[{}..{}), ",
                    name, sw.first_sample_global, sw.last_sample_global
                ));
            }
            log::debug!(
                "⚠️  Sample {} not found in any subswath (beyond tolerance). Available ranges: {}",
                sample,
                subswath_ranges
            );
        }

        None
    }

    /// Get the slant_range_time for a given native range sample
    ///
    /// Uses subswath-specific slant_range_time when available, falls back to
    /// the global params.slant_range_time when subswath cannot be determined.
    #[inline]
    pub fn get_slant_range_time_for_sample(&self, native_range_sample: f64) -> f64 {
        if let Some(subswath) = self.find_subswath_for_sample(native_range_sample) {
            subswath.slant_range_time
        } else {
            // Fallback to global slant_range_time
            self.slant_range_time
        }
    }

    /// Compute range pixel from slant range using per-subswath geometry
    ///
    /// This is the key function that enables correct geocoding for merged IW data
    /// by using the appropriate slant_range_time for each subswath.
    ///
    /// For merged IW data, each subswath has its own slant_range_time which is the
    /// two-way travel time to sample 0 of that subswath's local coordinate system.
    /// The global sample index = first_sample_global + local_sample.
    #[inline]
    pub fn slant_range_to_native_pixel(&self, slant_range: f64) -> f64 {
        let two_way_time = 2.0 * slant_range / self.speed_of_light;
        let range_pixel_spacing_time = 2.0 * self.range_pixel_spacing / self.speed_of_light;

        // For IW mode with subswaths, we need to find which subswath this slant range falls into
        // Each subswath covers a different range of slant range times
        if !self.subswaths.is_empty() {
            // Find subswath by slant range time, not by pixel (which we don't know yet)
            // Sort subswaths by increasing first_sample_global so we prefer near-range (IW1) matches.
            let mut subs: Vec<&SubSwath> = self.subswaths.values().collect();
            subs.sort_by_key(|sw| sw.first_sample_global);

            // Store first subswath for fallback case
            let first_subswath = subs.first();

            // IMPROVED tolerance for subswath matching (5% of range pixel spacing time)
            // Increased from 1% to 5% to handle DEM elevation uncertainties and geoid conversion errors
            // DEM errors can cause 50-100m slant range errors, which translates to ~0.3-0.7ms timing error
            // 5% tolerance = ~0.15ms, which is reasonable for edge cases
            let tolerance = range_pixel_spacing_time * 0.05;

            // Track best match for fallback (nearest subswath by slant_range_time)
            let mut best_match: Option<(&SubSwath, f64)> = None;
            let mut best_distance = f64::INFINITY;

            for subswath in &subs {
                let subswath_srt = subswath.slant_range_time;
                let subswath_num_samples = subswath.range_samples;

                // Slant range time at the end of this subswath
                let subswath_end_srt =
                    subswath_srt + (subswath_num_samples as f64) * range_pixel_spacing_time;

                // Check if two_way_time falls within this subswath's range (with tolerance)
                if two_way_time >= (subswath_srt - tolerance)
                    && two_way_time <= (subswath_end_srt + tolerance)
                {
                    // Local pixel within this subswath
                    let local_pixel = (two_way_time - subswath_srt) / range_pixel_spacing_time;
                    // Global pixel = first_sample_global + local_pixel
                    let global_pixel = subswath.first_sample_global as f64 + local_pixel;
                    return global_pixel;
                }

                // Track distance to this subswath for fallback
                let distance = if two_way_time < subswath_srt {
                    subswath_srt - two_way_time
                } else if two_way_time > subswath_end_srt {
                    two_way_time - subswath_end_srt
                } else {
                    0.0 // Should have matched above, but handle edge case
                };

                if distance < best_distance {
                    best_distance = distance;
                    best_match = Some((subswath, distance));
                }
            }

            // IMPROVED FALLBACK LOGIC: If no exact match, use best match if within reasonable tolerance
            // This handles cases where computed slant ranges are slightly off due to:
            // - DEM elevation uncertainties (can cause 50-100m errors in slant range)
            // - Geoid conversion errors (2-5m vertical error → ~10-25m slant range error)
            // - Coordinate transformation edge cases
            // - Numerical precision issues
            //
            // IMPROVED tolerance: 1.0ms = 150km slant range (at speed of light)
            // Increased from 0.5ms to 1.0ms to better handle geoid conversion errors and DEM uncertainties
            // Sentinel-1 IW subswaths are typically separated by ~20-25k range samples
            // (~46-58km slant range), so 150km tolerance is larger than subswath spacing
            // but necessary to handle combined DEM + geoid errors and edge cases.
            // This is safe because we still prefer exact matches when available.
            const MAX_FALLBACK_TOLERANCE: f64 = 0.001; // 1.0ms tolerance (~150km slant range)

            if let Some((best_subswath, distance)) = best_match {
                if distance <= MAX_FALLBACK_TOLERANCE {
                    // Use best match subswath
                    let local_pixel =
                        (two_way_time - best_subswath.slant_range_time) / range_pixel_spacing_time;
                    let global_pixel = best_subswath.first_sample_global as f64 + local_pixel;

                    // Log diagnostic information (only first few times to avoid spam)
                    static FALLBACK_LOG_COUNT: std::sync::atomic::AtomicUsize =
                        std::sync::atomic::AtomicUsize::new(0);
                    let count =
                        FALLBACK_LOG_COUNT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if count < 5 {
                        log::debug!(
                            "⚠️  Subswath fallback match: two_way_time={:.9}s, best_subswath_srt={:.9}s, distance={:.9}s, global_pixel={:.1}",
                            two_way_time, best_subswath.slant_range_time, distance, global_pixel
                        );
                    }

                    return global_pixel;
                }
            }

            // IMPROVED FALLBACK: Handle cases where two_way_time is less than subswath slant_range_time
            // This can occur when:
            // - DEM elevation is too high (reducing computed slant range)
            // - Edge pixels fall slightly outside subswath bounds
            // - Numerical precision issues
            //
            // Strategy: Use the first subswath (IW1, near-range) if two_way_time is less than
            // any subswath's slant_range_time, but clamp the result to prevent negative pixels.
            if let Some(first_subswath) = first_subswath {
                // Check if two_way_time is less than the first subswath's slant_range_time
                // or less than the global slant_range_time (which is typically from IW1)
                let min_srt = first_subswath.slant_range_time.min(self.slant_range_time);
                if two_way_time < min_srt {
                    // Calculate local pixel (may be negative)
                    let local_pixel =
                        (two_way_time - first_subswath.slant_range_time) / range_pixel_spacing_time;

                    // Clamp local_pixel to prevent negative global pixels
                    // Allow small negative values (within 10 pixels) for edge cases
                    let clamped_local_pixel = local_pixel.max(-10.0);
                    let global_pixel =
                        first_subswath.first_sample_global as f64 + clamped_local_pixel;

                    // Ensure global_pixel is non-negative
                    let final_global_pixel = global_pixel.max(0.0);

                    // Log diagnostic information
                    static FIRST_SUBSWATH_FALLBACK_COUNT: std::sync::atomic::AtomicUsize =
                        std::sync::atomic::AtomicUsize::new(0);
                    let count = FIRST_SUBSWATH_FALLBACK_COUNT
                        .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if count < 10 {
                        if local_pixel < 0.0 {
                            log::warn!(
                                "⚠️  First subswath fallback (clamped): two_way_time={:.9}s < min_srt={:.9}s, local_pixel={:.1} (clamped to {:.1}), global_pixel={:.1}",
                                two_way_time, min_srt, local_pixel, clamped_local_pixel, final_global_pixel
                            );
                        } else {
                            log::debug!(
                                "ℹ️  First subswath fallback: two_way_time={:.9}s < min_srt={:.9}s, global_pixel={:.1}",
                                two_way_time, min_srt, final_global_pixel
                            );
                        }
                    }

                    return final_global_pixel;
                }

                // Also handle case where two_way_time is slightly less than first subswath's slant_range_time
                // but greater than global slant_range_time (edge case)
                let tolerance_srt = range_pixel_spacing_time * 10.0; // 10 pixel tolerance
                if two_way_time >= (first_subswath.slant_range_time - tolerance_srt)
                    && two_way_time < first_subswath.slant_range_time
                {
                    // Use first subswath with small negative local pixel (clamped)
                    let local_pixel =
                        (two_way_time - first_subswath.slant_range_time) / range_pixel_spacing_time;
                    let clamped_local_pixel = local_pixel.max(-5.0); // Allow up to 5 pixels negative
                    let global_pixel =
                        (first_subswath.first_sample_global as f64 + clamped_local_pixel).max(0.0);

                    static EDGE_FALLBACK_COUNT: std::sync::atomic::AtomicUsize =
                        std::sync::atomic::AtomicUsize::new(0);
                    let count =
                        EDGE_FALLBACK_COUNT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if count < 10 {
                        log::debug!(
                            "ℹ️  Edge case fallback: two_way_time={:.9}s near first_subswath_srt={:.9}s, global_pixel={:.1}",
                            two_way_time, first_subswath.slant_range_time, global_pixel
                        );
                    }

                    return global_pixel;
                }
            }

            // If we didn't find a matching subswath and two_way_time >= global slant_range_time,
            // fall through to global calculation. This can happen in overlap regions or at edges.
        }

        // Fallback: use global slant_range_time
        (two_way_time - self.slant_range_time) / range_pixel_spacing_time
    }

    /// Validate if a native range pixel is within valid bounds
    ///
    /// Checks against both global bounds and subswath-specific valid sample ranges.
    #[inline]
    pub fn is_valid_range_pixel(&self, native_range_pixel: f64) -> bool {
        if native_range_pixel < 0.0 {
            return false;
        }

        // Check against subswath valid sample ranges if available
        if let Some(subswath) = self.find_subswath_for_sample(native_range_pixel) {
            if let (Some(first_valid), Some(last_valid)) =
                (subswath.valid_first_sample, subswath.valid_last_sample)
            {
                let sample = native_range_pixel as usize;
                return sample >= first_valid && sample <= last_valid;
            }
        }

        // Fallback: check against global first/last valid sample
        if let (Some(first), Some(last)) = (self.first_valid_sample, self.last_valid_sample) {
            let sample = native_range_pixel as usize;
            return sample >= first && sample <= last;
        }

        // No explicit bounds, accept any non-negative pixel
        true
    }
}
