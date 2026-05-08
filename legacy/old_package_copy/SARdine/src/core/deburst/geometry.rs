#![allow(dead_code)]
use crate::constants::physical::SPEED_OF_LIGHT_M_S;
use crate::io::annotation::raw::SwathTiming;
use crate::types::{SarError, SarResult, SubSwath};
use chrono::NaiveDateTime;

/// Holds per-line timing derived from annotation (preferred over spacing/velocity)
#[derive(Clone, Copy)]
pub(crate) struct LineTiming {
    /// seconds relative to polynomial reference time (annotation-aware)
    pub(crate) t_az: f64,
}

/// Range-dependent polynomial model sampled at discrete slant-range times
#[derive(Debug, Clone)]
pub(crate) struct RangePolyEntry {
    pub(crate) range_time: f64,
    pub(crate) coeffs: Vec<f64>,
    pub(crate) t0: Option<f64>,
}

/// Range-dependent polynomial model sampled at discrete slant-range times
#[derive(Debug, Clone)]
pub(crate) struct RangePolynomial {
    pub(crate) entries: Vec<RangePolyEntry>,
}

// REMOVED: build_line_timing (use build_line_timing_with_offset for proper reference time handling)

/// Build per-line timings with correct annotation reference time
pub(crate) fn build_line_timing_with_offset(
    lines: usize,
    az_time_interval_s: f64,
    time_offset_s: f64,
) -> Vec<LineTiming> {
    let center = (lines as f64 - 1.0) * 0.5;
    (0..lines)
        .map(|l| {
            let t_burst_relative = (l as f64 - center) * az_time_interval_s;
            LineTiming {
                t_az: t_burst_relative + time_offset_s,
            }
        })
        .collect()
}

/// Evaluate Doppler centroid (Hz) and FM rate (Hz/s) using slant range time τ with reference t0
pub(crate) fn eval_dc_fm(
    range_time_s: f64,
    dc_poly: &[f64],
    fm_poly: &[f64],
    dc_t0: Option<f64>,
    fm_t0: Option<f64>,
) -> (f64, f64) {
    let p = |c: &[f64], x: f64| c.iter().rev().fold(0.0, |acc, &a| acc * x + a);
    let tau_dc = range_time_s - dc_t0.unwrap_or(0.0);
    let tau_fm = range_time_s - fm_t0.unwrap_or(0.0);
    (p(dc_poly, tau_dc), p(fm_poly, tau_fm))
}

/// Evaluate Doppler centroid and FM rate with RANGE DEPENDENCY for TOPS IW alignment
pub(crate) fn eval_dc_fm_2d(
    range_sample: usize,
    dc_poly: &[f64],
    fm_poly: &[f64],
    dc_range_poly: Option<&RangePolynomial>,
    fm_range_poly: Option<&RangePolynomial>,
    range_pixel_spacing: f64,
    slant_range_time: f64,
    dc_t0: Option<f64>,
    fm_t0: Option<f64>,
) -> (f64, f64) {
    let range_time =
        slant_range_time + (range_sample as f64 * range_pixel_spacing * 2.0 / SPEED_OF_LIGHT_M_S);
    let tau_dc = range_time - dc_t0.unwrap_or(0.0);
    let tau_fm = range_time - fm_t0.unwrap_or(0.0);

    let f_dc = if let Some(model) = dc_range_poly {
        model.evaluate(range_time, dc_t0)
    } else if dc_poly.len() > 3 {
        eval_poly_2d(dc_poly, tau_dc, range_time)
    } else {
        let p = |c: &[f64], x: f64| c.iter().rev().fold(0.0, |acc, &a| acc * x + a);
        p(dc_poly, tau_dc)
    };

    let fm_rate = if let Some(model) = fm_range_poly {
        model.evaluate(range_time, fm_t0)
    } else if fm_poly.len() > 3 {
        eval_poly_2d(fm_poly, tau_fm, range_time)
    } else {
        let p = |c: &[f64], x: f64| c.iter().rev().fold(0.0, |acc, &a| acc * x + a);
        p(fm_poly, tau_fm)
    };

    (f_dc, fm_rate)
}

/// Evaluate 2D polynomial: f(t, r) = sum(c_ij × t^i × r^j)
pub(crate) fn eval_poly_2d(coeffs: &[f64], t: f64, r: f64) -> f64 {
    if coeffs.len() <= 3 {
        return coeffs.iter().rev().fold(0.0, |acc, &c| acc * t + c);
    }

    let mut result = 0.0;
    let max_order = ((coeffs.len() as f64).sqrt().ceil() as usize).max(1);

    for (idx, &coeff) in coeffs.iter().enumerate() {
        let i = idx / max_order;
        let j = idx % max_order;
        result += coeff * t.powi(i as i32) * r.powi(j as i32);
    }
    result
}

impl RangePolynomial {
    pub(crate) fn from_dc_estimates(estimates: &[DcEstimate]) -> Option<Self> {
        Self::from_estimates(estimates, |e| e.slant_range_time, |e| e.t0, |e| &e.coeffs)
    }

    pub(crate) fn from_fm_estimates(estimates: &[FmEstimate]) -> Option<Self> {
        Self::from_estimates(estimates, |e| e.slant_range_time, |e| e.t0, |e| &e.coeffs)
    }

    fn from_estimates<T, FRange, FT0, FCoeff>(
        estimates: &[T],
        range_fn: FRange,
        t0_fn: FT0,
        coeff_fn: FCoeff,
    ) -> Option<Self>
    where
        FRange: Fn(&T) -> Option<f64>,
        FT0: Fn(&T) -> Option<f64>,
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
                    t0: t0_fn(e),
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

    pub(crate) fn evaluate(&self, range_time: f64, default_t0: Option<f64>) -> f64 {
        if self.entries.is_empty() {
            return 0.0;
        }

        let eval_entry = |entry: &RangePolyEntry| {
            let t0 = entry.t0.or(default_t0).unwrap_or(0.0);
            let tau = range_time - t0;
            eval_poly_2d(&entry.coeffs, tau, range_time)
        };

        if self.entries.len() == 1 {
            return eval_entry(&self.entries[0]);
        }

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

    pub(crate) fn samples(&self) -> usize {
        self.entries.len()
    }
}

/// TOPSAR Burst information for proper debursting
#[derive(Debug, Clone)]
pub struct BurstInfo {
    pub burst_id: usize,
    pub start_line: usize,
    pub end_line: usize,
    pub start_sample: usize,
    pub end_sample: usize,
    pub azimuth_time: String,
    pub sensing_time: String,
    pub first_valid_sample: Vec<i32>,
    pub last_valid_sample: Vec<i32>,
    pub byte_offset: u64,
    /// Source subswath identifier (e.g., IW1/IW2/IW3) used to tripwire metadata mixups
    pub source_subswath: Option<String>,
    pub azimuth_fm_rate: f64,
    pub azimuth_steering_rate: f64,
    pub slant_range_time: f64,
    pub doppler_centroid: f64,
    pub azimuth_bandwidth: f64,
    pub range_sampling_rate: f64,
    pub range_pixel_spacing: f64,
    pub azimuth_pixel_spacing: f64,
    pub azimuth_time_interval: f64,
    pub dc_polynomial: Vec<f64>,
    pub fm_polynomial: Vec<f64>,
    pub(crate) dc_range_poly: Option<RangePolynomial>,
    pub(crate) fm_range_poly: Option<RangePolynomial>,
    pub dc_polynomial_t0: Option<f64>,
    pub fm_polynomial_t0: Option<f64>,
    pub burst_reference_time_seconds: Option<f64>,
}

impl BurstInfo {
    pub fn lines(&self) -> usize {
        self.end_line.saturating_sub(self.start_line) + 1
    }

    fn eval_dc_poly(poly: &[f64], t_az: f64) -> f64 {
        poly.iter()
            .enumerate()
            .fold(0.0, |acc, (i, &coeff)| acc + coeff * t_az.powi(i as i32))
    }

    pub fn valid_samples_for_line(&self, line: usize) -> (usize, usize) {
        if line >= self.first_valid_sample.len() || line >= self.last_valid_sample.len() {
            return (0, 0);
        }

        let first = self.first_valid_sample[line].max(0) as usize;
        let last = self.last_valid_sample[line].max(0) as usize;

        if first <= last {
            (first, last)
        } else {
            (0, 0)
        }
    }

    pub fn with_enhanced_timing(
        mut self,
        azimuth_time_interval: f64,
        dc_polynomial: Vec<f64>,
        fm_polynomial: Vec<f64>,
    ) -> Self {
        self.azimuth_time_interval = azimuth_time_interval;
        self.dc_polynomial = dc_polynomial;
        self.fm_polynomial = fm_polynomial;
        self
    }
}

/// DC estimate parsed from annotation (coefficients + optional reference time)
#[derive(Debug, Clone)]
pub(crate) struct DcEstimate {
    pub(crate) coeffs: Vec<f64>,
    pub(crate) t0: Option<f64>,
    pub(crate) azimuth_time: Option<f64>,
    pub(crate) slant_range_time: Option<f64>,
}

/// FM estimate parsed from annotation (coefficients + optional reference time)
#[derive(Debug, Clone)]
pub(crate) struct FmEstimate {
    pub(crate) coeffs: Vec<f64>,
    pub(crate) t0: Option<f64>,
    pub(crate) azimuth_time: Option<f64>,
    pub(crate) slant_range_time: Option<f64>,
}

/// Derived burst geometry (moved from io::annotation::derive)
#[derive(Debug, Clone)]
pub(crate) struct DerivedBurstGeometry {
    pub(crate) total_bursts: usize,
    pub(crate) valid_bursts: usize,
    pub(crate) lines_per_burst: Option<usize>,
    pub(crate) samples_per_burst: Option<usize>,
    pub(crate) total_lines: usize,
    pub(crate) first_valid_line: Option<usize>,
    pub(crate) last_valid_line_exclusive: Option<usize>,
    pub(crate) first_valid_sample: Option<usize>,
    pub(crate) last_valid_sample_exclusive: Option<usize>,
}

pub(crate) fn derive_burst_geometry(swath_timing: &SwathTiming) -> Option<DerivedBurstGeometry> {
    let bursts = swath_timing
        .burst_list
        .as_ref()
        .and_then(|list| list.bursts.as_ref())?;

    if bursts.is_empty() {
        return None;
    }

    let mut line_cursor = 0usize;
    let mut first_valid_line = None;
    let mut last_valid_line = None;
    let mut first_valid_sample = None;
    let mut last_valid_sample = None;
    let mut valid_bursts = 0usize;

    for burst in bursts {
        let lines_in_burst = burst
            .first_valid_sample
            .len()
            .max(burst.last_valid_sample.len());
        let mut burst_has_valid = false;

        for line_idx in 0..lines_in_burst {
            let first_sample = burst
                .first_valid_sample
                .get(line_idx)
                .copied()
                .unwrap_or(-1);
            let last_sample = burst.last_valid_sample.get(line_idx).copied().unwrap_or(-1);

            if first_sample >= 0 && last_sample >= first_sample {
                let global_line = line_cursor + line_idx;
                first_valid_line =
                    Some(first_valid_line.map_or(global_line, |v: usize| v.min(global_line)));
                last_valid_line = Some(
                    last_valid_line.map_or(global_line + 1, |v: usize| v.max(global_line + 1)),
                );

                let first_usize = first_sample as usize;
                let last_usize = (last_sample as usize).saturating_add(1);
                first_valid_sample =
                    Some(first_valid_sample.map_or(first_usize, |v: usize| v.min(first_usize)));
                last_valid_sample =
                    Some(last_valid_sample.map_or(last_usize, |v: usize| v.max(last_usize)));

                burst_has_valid = true;
            }
        }

        if burst_has_valid {
            valid_bursts += 1;
        }

        line_cursor += lines_in_burst;
    }

    Some(DerivedBurstGeometry {
        total_bursts: bursts.len(),
        valid_bursts,
        lines_per_burst: swath_timing.lines_per_burst.map(|v| v as usize),
        samples_per_burst: swath_timing.samples_per_burst.map(|v| v as usize),
        total_lines: line_cursor,
        first_valid_line,
        last_valid_line_exclusive: last_valid_line,
        first_valid_sample,
        last_valid_sample_exclusive: last_valid_sample,
    })
}

pub fn extract_burst_info_from_annotation(
    annotation_data: &str,
    total_lines: usize,
    total_samples: usize,
) -> SarResult<Vec<BurstInfo>> {
    extract_burst_info_from_annotation_with_subswath(
        annotation_data,
        total_lines,
        total_samples,
        None,
    )
}

pub fn extract_burst_info_from_annotation_with_subswath(
    annotation_data: &str,
    total_lines: usize,
    total_samples: usize,
    subswath: Option<&SubSwath>,
) -> SarResult<Vec<BurstInfo>> {
    log::info!("🔍 Extracting burst information from Sentinel-1 annotation");
    log::debug!(
        "Total image dimensions: {} x {}",
        total_lines,
        total_samples
    );

    let base_line = subswath
        .map(|sw| sw.first_line_global as f64)
        .unwrap_or(0.0)
        .max(0.0) as usize;

    // CRITICAL FIX: Burst sample positions must stay in LOCAL subswath coordinates.
    // The SubSwath.valid_first_sample was incorrectly set to the GLOBAL sample offset
    // (e.g., 19934 for IW2), causing severe width truncation in deburst when combined
    // with local total_samples (25957 - 19934 = only 6023 samples).
    //
    // Global offsets are ONLY relevant for merge (aligning subswaths in the output grid).
    // During deburst, bursts are indexed in local coordinates (0..samples_per_burst).
    // Force base_sample = 0 to preserve full subswath width.
    let base_sample = 0usize;

    match parse_topsar_burst_info_with_subswath(
        annotation_data,
        total_lines,
        total_samples,
        subswath,
    ) {
        Ok(bursts) if !bursts.is_empty() => {
            log::info!(
                "✅ Successfully extracted {} TOPSAR bursts from annotation",
                bursts.len()
            );

            let mut adjusted: Vec<BurstInfo> = Vec::with_capacity(bursts.len());

            for (idx, mut burst) in bursts.into_iter().enumerate() {
                let burst_lines = burst.first_valid_sample.len();
                if burst_lines == 0 {
                    return Err(SarError::Processing(format!(
                        "Burst {} has no first_valid_sample entries",
                        idx
                    )));
                }

                let rel_start_line = idx * burst_lines;
                burst.start_line = base_line + rel_start_line;

                if burst.start_line >= total_lines {
                    return Err(SarError::Processing(format!(
                        "Burst {} start {} exceeds SLC height {}",
                        idx, burst.start_line, total_lines
                    )));
                }

                let lines_for_burst = burst_lines.min(total_lines - burst.start_line);
                if lines_for_burst == 0 {
                    return Err(SarError::Processing(format!(
                        "Burst {} has no remaining lines within SLC height {}",
                        idx, total_lines
                    )));
                }
                if lines_for_burst < burst_lines {
                    log::warn!(
                        "Burst {} truncated from {} to {} lines to fit SLC height {}",
                        idx,
                        burst_lines,
                        lines_for_burst,
                        total_lines
                    );
                }

                burst.end_line = burst
                    .start_line
                    .saturating_add(lines_for_burst.saturating_sub(1));

                let local_width = burst.end_sample.saturating_sub(burst.start_sample) + 1;
                let max_width = total_samples.saturating_sub(base_sample);
                let width = local_width.min(max_width.max(1));
                burst.start_sample = base_sample;
                burst.end_sample = base_sample + width.saturating_sub(1);

                let width = burst.end_sample.saturating_sub(burst.start_sample) + 1;

                let clamp_sample = |v: i32| -> i32 {
                    let mut adj = v - base_sample as i32;
                    if adj < 0 {
                        adj = 0;
                    }
                    let max_val = (width as i32).saturating_sub(1);
                    if adj > max_val {
                        adj = max_val;
                    }
                    adj
                };

                if burst.first_valid_sample.len() != lines_for_burst {
                    let fill = *burst.first_valid_sample.last().unwrap_or(&0);
                    burst
                        .first_valid_sample
                        .truncate(lines_for_burst.min(burst.first_valid_sample.len()));
                    burst.first_valid_sample.resize(lines_for_burst, fill);
                }

                if burst.last_valid_sample.len() != lines_for_burst {
                    let fill = *burst.last_valid_sample.last().unwrap_or(&0);
                    burst
                        .last_valid_sample
                        .truncate(lines_for_burst.min(burst.last_valid_sample.len()));
                    burst.last_valid_sample.resize(lines_for_burst, fill);
                }

                burst.first_valid_sample = burst
                    .first_valid_sample
                    .iter()
                    .map(|&v| clamp_sample(v))
                    .collect();

                burst.last_valid_sample = burst
                    .last_valid_sample
                    .iter()
                    .map(|&v| clamp_sample(v))
                    .collect();

                adjusted.push(burst);
            }

            Ok(adjusted)
        }
        Ok(_) => Err(SarError::Processing(
            "No valid bursts found in annotation data".to_string(),
        )),
        Err(e) => Err(SarError::Processing(format!(
            "Failed to parse burst information: {}",
            e
        ))),
    }
}

pub(crate) fn parse_topsar_burst_info(
    annotation_data: &str,
    total_lines: usize,
    total_samples: usize,
) -> SarResult<Vec<BurstInfo>> {
    parse_topsar_burst_info_with_subswath(annotation_data, total_lines, total_samples, None)
}

pub(crate) fn parse_topsar_burst_info_with_subswath(
    annotation_data: &str,
    total_lines: usize,
    total_samples: usize,
    subswath: Option<&SubSwath>,
) -> SarResult<Vec<BurstInfo>> {
    log::info!("🎯 Parsing TOPSAR burst information with enhanced parameters");

    let mut burst_info = Vec::new();

    let azimuth_fm_rate = extract_parameter_string(annotation_data, "<azimuthFmRatePolynomial", "</azimuthFmRatePolynomial>")
        .and_then(|s| {
            if let Some(content_start) = s.find('>') {
                let content = &s[content_start + 1..];
                let coeffs: Vec<&str> = content.split_whitespace().collect();
                coeffs.first().and_then(|s| s.parse::<f64>().ok())
            } else {
                let coeffs: Vec<&str> = s.split_whitespace().collect();
                coeffs.first().and_then(|s| s.parse::<f64>().ok())
            }
        })
        .or_else(|| extract_parameter(annotation_data, "<azimuthFmRate>", "</azimuthFmRate>"))
        .or_else(|| extract_parameter(annotation_data, "<azimuthFMRate>", "</azimuthFMRate>"))
        .ok_or_else(|| SarError::Processing("❌ SCIENTIFIC ERROR: Azimuth FM rate not found in annotation XML! Real Sentinel-1 parameters required - no fallbacks allowed.".to_string()))?;
    let azimuth_steering_rate = extract_parameter(annotation_data, "<azimuthSteeringRate>", "</azimuthSteeringRate>")
        .ok_or_else(|| SarError::Processing("❌ SCIENTIFIC ERROR: Azimuth steering rate not found in annotation XML! Real Sentinel-1 parameters required - no fallbacks allowed.".to_string()))?;
    let range_sampling_rate = extract_parameter(annotation_data, "<rangeSamplingRate>", "</rangeSamplingRate>")
        .ok_or_else(|| SarError::Processing("❌ SCIENTIFIC ERROR: Range sampling rate not found in annotation XML! Real Sentinel-1 parameters required - no fallbacks allowed.".to_string()))?;

    let range_pixel_spacing = extract_parameter(annotation_data, "<rangePixelSpacing>", "</rangePixelSpacing>")
        .ok_or_else(|| SarError::Processing("❌ SCIENTIFIC ERROR: Range pixel spacing not found in annotation XML! Real Sentinel-1 annotation required - no synthetic fallbacks allowed for research-grade processing.".to_string()))?;
    let azimuth_pixel_spacing = extract_parameter(annotation_data, "<azimuthPixelSpacing>", "</azimuthPixelSpacing>")
        .ok_or_else(|| SarError::Processing("❌ SCIENTIFIC ERROR: Azimuth pixel spacing not found in annotation XML! Real Sentinel-1 annotation required - no synthetic fallbacks allowed for research-grade processing.".to_string()))?;

    let slant_range_time = extract_parameter(annotation_data, "<slantRangeTime>", "</slantRangeTime>")
        .or_else(|| subswath.map(|sw| sw.slant_range_time))
        .ok_or_else(|| SarError::Processing("❌ SCIENTIFIC ERROR: Slant range time not found in annotation XML or subswath metadata.".to_string()))?;

    let lines_per_burst = extract_parameter(annotation_data, "<linesPerBurst>", "</linesPerBurst>")
        .ok_or_else(|| SarError::Processing("❌ SCIENTIFIC ERROR: Lines per burst not found in annotation XML! Real Sentinel-1 burst parameters required - no synthetic values allowed.".to_string()))? as usize;

    let samples_per_burst = extract_parameter(annotation_data, "<samplesPerBurst>", "</samplesPerBurst>")
        .ok_or_else(|| SarError::Processing("❌ SCIENTIFIC ERROR: Samples per burst not found in annotation XML! Real Sentinel-1 burst parameters required - no synthetic values allowed.".to_string()))? as usize;

    let azimuth_time_interval = extract_parameter(
        annotation_data,
        "<azimuthTimeInterval>",
        "</azimuthTimeInterval>",
    )
    .ok_or_else(|| {
        SarError::Processing(
            "CRITICAL: azimuthTimeInterval missing from annotation; DC-aware deburst requires annotation-derived azimuth timing (no PRF fallback)."
                .to_string(),
        )
    })?;

    let dc_estimates = extract_dc_estimates_from_annotation(annotation_data);

    let default_dc_poly = if let Some(sw) = subswath {
        if let Some(ref dc_poly) = sw.dc_polynomial {
            log::info!(
                "✅ Using pre-parsed DC polynomial from SubSwath with {} coefficients: {:?}",
                dc_poly.len(),
                dc_poly
            );
            Some(dc_poly.clone())
        } else {
            log::warn!(
                "⚠️ SubSwath provided but dc_polynomial is None, will attempt XML-derived DC polynomials"
            );
            None
        }
    } else {
        log::warn!("⚠️ No SubSwath provided, will parse DC polynomials from annotation");
        None
    };

    let fm_estimates = extract_fm_estimates_from_annotation(annotation_data);

    let dc_range_polynomial: Option<RangePolynomial> = dc_estimates
        .as_ref()
        .and_then(|v| RangePolynomial::from_dc_estimates(v.as_slice()));
    if let Some(model) = &dc_range_polynomial {
        log::info!(
            "✅ Detected range-dependent DC polynomial grid ({} samples)",
            model.samples()
        );
    }

    let fm_range_polynomial: Option<RangePolynomial> = fm_estimates
        .as_ref()
        .and_then(|v| RangePolynomial::from_fm_estimates(v.as_slice()));
    if let Some(model) = &fm_range_polynomial {
        log::info!(
            "✅ Detected range-dependent FM polynomial grid ({} samples)",
            model.samples()
        );
    }

    let default_fm_poly = extract_polynomial_coefficients(
        annotation_data,
        "<fmRatePolynomial>",
        "</fmRatePolynomial>",
    )
    .or_else(|| {
        log::info!("Using azimuth FM rate as constant polynomial");
        Some(vec![azimuth_fm_rate, 0.0])
    });

    log::info!(
        "📊 TOPSAR parameters: lines_per_burst={}, samples_per_burst={}",
        lines_per_burst,
        samples_per_burst
    );
    log::info!("📊 Range sampling rate: {:.0} Hz", range_sampling_rate);
    log::info!(
        "📊 Azimuth steering rate: {:.6} rad/s",
        azimuth_steering_rate
    );

    let burst_pattern = regex::Regex::new(
        r"(?s)<burst>.*?<azimuthTime>([^<]+)</azimuthTime>.*?<sensingTime>([^<]+)</sensingTime>.*?<byteOffset>([^<]+)</byteOffset>.*?<firstValidSample[^>]*>([^<]+)</firstValidSample>.*?<lastValidSample[^>]*>([^<]+)</lastValidSample>.*?</burst>"
    ).map_err(|e| SarError::Processing(format!("Failed to compile burst regex: {}", e)))?;

    let burst_matches: Vec<_> = burst_pattern.captures_iter(annotation_data).collect();

    if burst_matches.is_empty() {
        log::error!("❌ No burst information found with regex pattern");

        if let Ok(count_regex) = regex::Regex::new(r#"<burstList count="(\d+)">"#) {
            if let Some(count_match) = count_regex.captures(annotation_data) {
                if let Ok(burst_count) = count_match[1].parse::<usize>() {
                    log::info!(
                        "📊 Found burstList with {} bursts, but couldn't parse individual bursts",
                        burst_count
                    );
                }
            }
        }

        return Err(SarError::Processing(
            "No burst information found in annotation".to_string(),
        ));
    }

    log::info!(
        "✅ Found {} burst matches in annotation",
        burst_matches.len()
    );

    for (i, captures) in burst_matches.iter().enumerate() {
        let azimuth_time = captures
            .get(1)
            .ok_or_else(|| {
                SarError::Processing(format!("Missing azimuth_time in burst {} regex match", i))
            })?
            .as_str()
            .to_string();
        let sensing_time = captures
            .get(2)
            .ok_or_else(|| {
                SarError::Processing(format!("Missing sensing_time in burst {} regex match", i))
            })?
            .as_str()
            .to_string();
        let byte_offset = captures
            .get(3)
            .ok_or_else(|| {
                SarError::Processing(format!("Missing byte_offset in burst {} regex match", i))
            })?
            .as_str()
            .parse::<u64>()
            .map_err(|e| {
                SarError::Processing(format!(
                    "Failed to parse byte_offset for burst {}: {}",
                    i, e
                ))
            })?;

        let first_valid_sample = parse_sample_array(
            captures
                .get(4)
                .ok_or_else(|| {
                    SarError::Processing(format!(
                        "Missing first_valid_sample in burst {} regex match",
                        i
                    ))
                })?
                .as_str(),
        );
        let last_valid_sample = parse_sample_array(
            captures
                .get(5)
                .ok_or_else(|| {
                    SarError::Processing(format!(
                        "Missing last_valid_sample in burst {} regex match",
                        i
                    ))
                })?
                .as_str(),
        );

        if first_valid_sample.len() != lines_per_burst {
            log::warn!(
                "⚠️ Burst {}: firstValidSample length {} != lines_per_burst {}",
                i,
                first_valid_sample.len(),
                lines_per_burst
            );
        }

        let start_line = i * lines_per_burst;
        let end_line = ((i + 1) * lines_per_burst - 1).min(total_lines.saturating_sub(1));
        let start_sample = 0;
        let end_sample = samples_per_burst
            .saturating_sub(1)
            .min(total_samples.saturating_sub(1));

        log::info!(
            "📋 Burst {}: lines {}..{}, samples {}..{}, byte_offset={}",
            i,
            start_line,
            end_line,
            start_sample,
            end_sample,
            byte_offset
        );

        let burst_reference_time_seconds =
            parse_iso8601_seconds(&sensing_time).or_else(|| parse_iso8601_seconds(&azimuth_time));

        let dc_entry = dc_estimates.as_ref().and_then(|v| {
            let closest = v
                .iter()
                .filter_map(|e| e.slant_range_time.map(|rt| (rt, e)))
                .min_by(|a, b| {
                    (a.0 - slant_range_time)
                        .abs()
                        .partial_cmp(&(b.0 - slant_range_time).abs())
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(_, e)| e);

            closest.or_else(|| v.get(i)).or_else(|| v.first())
        });

        let fm_entry = fm_estimates.as_ref().and_then(|v| {
            let closest = v
                .iter()
                .filter_map(|e| e.slant_range_time.map(|rt| (rt, e)))
                .min_by(|a, b| {
                    (a.0 - slant_range_time)
                        .abs()
                        .partial_cmp(&(b.0 - slant_range_time).abs())
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(_, e)| e);

            closest.or_else(|| v.get(i)).or_else(|| v.first())
        });

        let dc_ref_epoch = dc_entry.and_then(|e| e.azimuth_time);
        let fm_ref_epoch = fm_entry.and_then(|e| e.azimuth_time);

        let dc_polynomial = if let Some(ref poly) = default_dc_poly {
            poly.clone()
        } else if let Some(entry) = dc_entry {
            entry.coeffs.clone()
        } else {
            return Err(SarError::Processing(
                "CRITICAL: DC polynomial not found in SubSwath or annotation. Expected <dopplerCentroid><dcEstimateList><dcEstimate><dataDcPolynomial>. ABORTING to prevent silent phase misalignment.".to_string(),
            ));
        };

        let mut dc_polynomial_t0 = dc_ref_epoch.or_else(|| dc_entry.and_then(|e| e.t0));
        if dc_polynomial_t0.map(|t0| t0.abs() < 1.0e6).unwrap_or(false) {
            if let Some(t_abs) = burst_reference_time_seconds {
                let adjusted = t_abs + dc_polynomial_t0.unwrap();
                log::warn!(
                    "⚠️  DC polynomial t0 looks relative ({:.6}); anchoring to burst absolute time → {:.6}",
                    dc_polynomial_t0.unwrap(),
                    adjusted
                );
                dc_polynomial_t0 = Some(adjusted);
            }
        }
        if dc_polynomial_t0.is_none() {
            if let Some(sw) = subswath {
                if let Some(t0) = sw.dc_polynomial_t0 {
                    log::info!("✅ Using SubSwath DC polynomial t0 from cached metadata");
                    dc_polynomial_t0 = Some(t0);
                }
            }

            if dc_polynomial_t0.is_none() {
                if let Some(t0) = extract_first_dc_t0(annotation_data) {
                    log::info!("✅ Extracted DC polynomial t0 directly from annotation (fallback)");
                    dc_polynomial_t0 = Some(t0);
                }
            }

            if dc_polynomial_t0.is_none() {
                if let Some(t0) = burst_reference_time_seconds {
                    log::warn!(
                        "⚠️  DC polynomial t0 missing from annotation + SubSwath; using burst sensing time ({:.6}s) as reference",
                        t0
                    );
                    dc_polynomial_t0 = Some(t0);
                }
            }
        }

        let fm_polynomial = if let Some(entry) = fm_entry {
            entry.coeffs.clone()
        } else if let Some(ref poly) = default_fm_poly {
            poly.clone()
        } else {
            vec![azimuth_fm_rate, 0.0]
        };

        let mut fm_polynomial_t0 = fm_ref_epoch
            .or_else(|| fm_entry.and_then(|e| e.t0))
            .or(dc_ref_epoch);
        if fm_polynomial_t0.map(|t0| t0.abs() < 1.0e6).unwrap_or(false) {
            if let Some(t_abs) = burst_reference_time_seconds {
                let adjusted = t_abs + fm_polynomial_t0.unwrap();
                log::warn!(
                    "⚠️  FM polynomial t0 looks relative ({:.6}); anchoring to burst absolute time → {:.6}",
                    fm_polynomial_t0.unwrap(),
                    adjusted
                );
                fm_polynomial_t0 = Some(adjusted);
            }
        }
        if fm_polynomial_t0.is_none() {
            if let Some(t0) = dc_polynomial_t0 {
                log::info!("✅ Using DC polynomial t0 as FM reference time");
                fm_polynomial_t0 = Some(t0);
            }

            if fm_polynomial_t0.is_none() {
                if let Some(t0) = burst_reference_time_seconds {
                    log::warn!(
                        "⚠️  FM polynomial t0 missing from annotation; using burst sensing time ({:.6}s) as reference",
                        t0
                    );
                    fm_polynomial_t0 = Some(t0);
                }
            }
        }

        if dc_polynomial_t0.is_none() || fm_polynomial_t0.is_none() {
            return Err(SarError::Processing(
                format!(
                    "CRITICAL: Missing DC/FM polynomial reference time (t0) for burst {} in annotation; no fallback allowed.",
                    i
                ),
            ));
        }

        burst_info.push(BurstInfo {
            burst_id: i,
            start_line,
            end_line,
            start_sample,
            end_sample,
            azimuth_time,
            sensing_time,
            first_valid_sample,
            last_valid_sample,
            byte_offset,
            azimuth_fm_rate,
            azimuth_steering_rate,
            slant_range_time,
            doppler_centroid: 0.0,
            azimuth_bandwidth: 320.0,
            range_sampling_rate,
            range_pixel_spacing,
            azimuth_pixel_spacing,
            azimuth_time_interval,
            dc_polynomial,
            fm_polynomial,
            dc_range_poly: dc_range_polynomial.clone(),
            fm_range_poly: fm_range_polynomial.clone(),
            dc_polynomial_t0,
            fm_polynomial_t0,
            burst_reference_time_seconds,
            source_subswath: subswath.map(|sw| sw.id.clone()),
        });
    }

    log::info!(
        "✅ Successfully parsed {} TOPSAR bursts with real parameters",
        burst_info.len()
    );
    Ok(burst_info)
}

pub(crate) fn parse_iso8601_seconds(value: &str) -> Option<f64> {
    let trimmed = value.trim();
    let fmt = "%Y-%m-%dT%H:%M:%S%.f";
    NaiveDateTime::parse_from_str(trimmed, fmt)
        .ok()
        .map(|dt| dt.and_utc())
        .map(|dt_utc| dt_utc.timestamp() as f64 + dt_utc.timestamp_subsec_nanos() as f64 * 1e-9)
}

fn extract_parameter(annotation_data: &str, start_tag: &str, end_tag: &str) -> Option<f64> {
    if let Some(start_pos) = annotation_data.find(start_tag) {
        let content_start = start_pos + start_tag.len();
        if let Some(end_pos) = annotation_data[content_start..].find(end_tag) {
            let content = &annotation_data[content_start..content_start + end_pos];
            return content.trim().parse::<f64>().ok();
        }
    }
    None
}

fn extract_parameter_string(
    annotation_data: &str,
    start_tag: &str,
    end_tag: &str,
) -> Option<String> {
    if let Some(start_pos) = annotation_data.find(start_tag) {
        let content_start = start_pos + start_tag.len();
        if let Some(end_pos) = annotation_data[content_start..].find(end_tag) {
            let content = &annotation_data[content_start..content_start + end_pos];
            return Some(content.trim().to_string());
        }
    }
    None
}

/// Extract DC estimates from annotation <dopplerCentroid><dcEstimateList>
pub(crate) fn extract_dc_estimates_from_annotation(
    annotation_data: &str,
) -> Option<Vec<DcEstimate>> {
    let block_re = regex::Regex::new(r"(?s)<dcEstimate[^>]*>(?P<body>.*?)</dcEstimate>").ok()?;
    let poly_re =
        regex::Regex::new(r"<dataDcPolynomial[^>]*>(?P<poly>[^<]+)</dataDcPolynomial>").ok()?;
    let t0_re = regex::Regex::new(r"<t0[^>]*>(?P<t0>[^<]+)</t0>").ok()?;
    let az_re = regex::Regex::new(r"<azimuthTime[^>]*>(?P<az>[^<]+)</azimuthTime>").ok()?;
    let srt_re = regex::Regex::new(r"<slantRangeTime[^>]*>(?P<srt>[^<]+)</slantRangeTime>").ok()?;

    let mut estimates: Vec<DcEstimate> = Vec::new();

    for caps in block_re.captures_iter(annotation_data) {
        let body = caps.name("body").map(|m| m.as_str()).unwrap_or("");

        let coeffs: Vec<f64> = poly_re
            .captures(body)
            .and_then(|m| m.name("poly"))
            .map(|m| m.as_str())
            .unwrap_or("")
            .split_whitespace()
            .filter_map(|s| s.parse::<f64>().ok())
            .collect();

        if coeffs.is_empty() {
            continue;
        }

        let t0 = t0_re
            .captures(body)
            .and_then(|m| m.name("t0"))
            .and_then(|m| m.as_str().trim().parse::<f64>().ok());

        let azimuth_time = az_re
            .captures(body)
            .and_then(|m| m.name("az"))
            .and_then(|m| parse_iso8601_seconds(m.as_str()));

        let slant_range_time = srt_re
            .captures(body)
            .and_then(|m| m.name("srt"))
            .and_then(|m| m.as_str().trim().parse::<f64>().ok());

        estimates.push(DcEstimate {
            coeffs,
            t0,
            azimuth_time,
            slant_range_time,
        });
    }

    if estimates.is_empty() {
        log::warn!("⚠️  dataDcPolynomial not found in annotation dcEstimateList");
        return None;
    }

    estimates.sort_by(|a, b| match (a.azimuth_time, b.azimuth_time) {
        (Some(at), Some(bt)) => at.partial_cmp(&bt).unwrap_or(std::cmp::Ordering::Equal),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => std::cmp::Ordering::Equal,
    });

    log::info!(
        "✅ Extracted {} DC estimates from annotation",
        estimates.len()
    );
    if let Some(first) = estimates.first() {
        log::info!(
            "   DC(t) = {} + {}*t + {}*t^2 + ... (t0={:?})",
            first.coeffs.get(0).unwrap_or(&0.0),
            first.coeffs.get(1).unwrap_or(&0.0),
            first.coeffs.get(2).unwrap_or(&0.0),
            first.t0
        );
        if estimates.iter().any(|e| e.slant_range_time.is_some()) {
            let mut ranges: Vec<f64> = estimates
                .iter()
                .filter_map(|e| e.slant_range_time)
                .collect();
            ranges.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            ranges.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
            log::info!("   Range-dependent DC grid with {} samples", ranges.len());
        }
    }

    Some(estimates)
}

/// Extract FM rate estimates from annotation <azimuthFmRateList><azimuthFmRate>
pub(crate) fn extract_fm_estimates_from_annotation(
    annotation_data: &str,
) -> Option<Vec<FmEstimate>> {
    let block_patterns = [
        r"(?s)<azimuthFmRate[^>]*>(?P<body>.*?)</azimuthFmRate>",
        r"(?s)<fmRate[^>]*>(?P<body>.*?)</fmRate>",
    ];

    let poly_patterns = [
        r"<azimuthFmRatePolynomial[^>]*>(?P<poly>[^<]+)</azimuthFmRatePolynomial>",
        r"<dataFmratePolynomial[^>]*>(?P<poly>[^<]+)</dataFmratePolynomial>",
    ];

    let t0_re = regex::Regex::new(r"<t0[^>]*>(?P<t0>[^<]+)</t0>").ok()?;
    let az_re = regex::Regex::new(r"<azimuthTime[^>]*>(?P<az>[^<]+)</azimuthTime>").ok()?;
    let srt_re = regex::Regex::new(r"<slantRangeTime[^>]*>(?P<srt>[^<]+)</slantRangeTime>").ok()?;

    let mut estimates: Vec<FmEstimate> = Vec::new();

    for block_pat in block_patterns {
        let block_re = match regex::Regex::new(block_pat) {
            Ok(r) => r,
            Err(_) => continue,
        };

        for caps in block_re.captures_iter(annotation_data) {
            let body = caps.name("body").map(|m| m.as_str()).unwrap_or("");

            let coeffs = poly_patterns
                .iter()
                .filter_map(|pat| regex::Regex::new(pat).ok())
                .find_map(|re| {
                    re.captures(body)
                        .and_then(|m| m.name("poly"))
                        .map(|m| m.as_str())
                })
                .map(|s| {
                    s.split_whitespace()
                        .filter_map(|v| v.parse::<f64>().ok())
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default();

            if coeffs.is_empty() {
                continue;
            }

            let t0 = t0_re
                .captures(body)
                .and_then(|m| m.name("t0"))
                .and_then(|m| m.as_str().trim().parse::<f64>().ok());

            let azimuth_time = az_re
                .captures(body)
                .and_then(|m| m.name("az"))
                .and_then(|m| parse_iso8601_seconds(m.as_str()));

            let slant_range_time = srt_re
                .captures(body)
                .and_then(|m| m.name("srt"))
                .and_then(|m| m.as_str().trim().parse::<f64>().ok());

            estimates.push(FmEstimate {
                coeffs,
                t0,
                azimuth_time,
                slant_range_time,
            });
        }
    }

    if estimates.is_empty() {
        log::warn!("⚠️  FM rate polynomial not found in annotation (azimuthFmRate or fmRate)");
        return None;
    }

    estimates.sort_by(|a, b| match (a.azimuth_time, b.azimuth_time) {
        (Some(at), Some(bt)) => at.partial_cmp(&bt).unwrap_or(std::cmp::Ordering::Equal),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => std::cmp::Ordering::Equal,
    });

    log::info!(
        "✅ Extracted {} FM estimates from annotation",
        estimates.len()
    );
    if let Some(first) = estimates.first() {
        log::info!(
            "   FM(t) = {} + {}*t + {}*t^2 + ... (t0={:?})",
            first.coeffs.get(0).unwrap_or(&0.0),
            first.coeffs.get(1).unwrap_or(&0.0),
            first.coeffs.get(2).unwrap_or(&0.0),
            first.t0
        );
        if estimates.iter().any(|e| e.slant_range_time.is_some()) {
            let mut ranges: Vec<f64> = estimates
                .iter()
                .filter_map(|e| e.slant_range_time)
                .collect();
            ranges.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            ranges.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
            log::info!("   Range-dependent FM grid with {} samples", ranges.len());
        }
    }

    Some(estimates)
}

fn extract_first_dc_t0(annotation_data: &str) -> Option<f64> {
    let re =
        regex::Regex::new(r"(?s)<dcEstimate[^>]*>.*?<t0[^>]*>(?P<t0>[^<]+)</t0>.*?</dcEstimate>")
            .ok()?;

    re.captures(annotation_data)
        .and_then(|caps| caps.name("t0"))
        .and_then(|m| m.as_str().trim().parse::<f64>().ok())
}

fn extract_polynomial_coefficients(
    annotation_data: &str,
    start_tag: &str,
    end_tag: &str,
) -> Option<Vec<f64>> {
    if let Some(content) = extract_parameter_string(annotation_data, start_tag, end_tag) {
        let coefficients: Vec<f64> = content
            .split_whitespace()
            .filter_map(|s| s.parse::<f64>().ok())
            .collect();
        if !coefficients.is_empty() {
            Some(coefficients)
        } else {
            None
        }
    } else {
        None
    }
}

fn parse_sample_array(data: &str) -> Vec<i32> {
    data.split_whitespace()
        .filter_map(|s| s.parse::<i32>().ok())
        .collect()
}
