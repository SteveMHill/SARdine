//! Quality metrics calculation for merge operations.

use super::super::overlap::OverlapRegion;
use super::super::types::QualityResults;
use crate::types::{SarError, SarResult};

/// Calculate quality metrics for the merge.
pub fn calculate_quality_metrics(overlaps: &[OverlapRegion]) -> SarResult<QualityResults> {
    if overlaps.is_empty() {
        return Err(SarError::Processing(
            "Cannot derive merge quality metrics without overlap regions".to_string(),
        ));
    }

    const PHASE_THRESHOLD: f32 = 0.85;
    const RADIOMETRIC_THRESHOLD: f32 = 0.90;
    const VALID_PIX_THRESHOLD: f32 = 0.90;
    const SEAM_THRESHOLD: f32 = 0.85;

    let mut warnings = Vec::new();
    let mut min_phase = 1.0_f32;
    let mut min_radiometric = 1.0_f32;
    let mut min_valid = 1.0_f32;
    let mut min_seam = 1.0_f32;

    for region in overlaps {
        let q = &region.quality_metrics;
        min_phase = min_phase.min(q.phase_coherence);
        min_radiometric = min_radiometric.min(q.radiometric_consistency);
        min_valid = min_valid.min(q.valid_pixel_percentage);
        min_seam = min_seam.min(q.seamline_quality);

        if q.phase_coherence < PHASE_THRESHOLD {
            warnings.push(format!(
                "Overlap {}-{} phase coherence {:.2} below {:.2}",
                region.swath1_id, region.swath2_id, q.phase_coherence, PHASE_THRESHOLD
            ));
        }
        if q.radiometric_consistency < RADIOMETRIC_THRESHOLD {
            warnings.push(format!(
                "Overlap {}-{} radiometric consistency {:.2} below {:.2}",
                region.swath1_id,
                region.swath2_id,
                q.radiometric_consistency,
                RADIOMETRIC_THRESHOLD
            ));
        }
        if q.valid_pixel_percentage < VALID_PIX_THRESHOLD {
            warnings.push(format!(
                "Overlap {}-{} valid coverage {:.1}% below {:.1}%",
                region.swath1_id,
                region.swath2_id,
                q.valid_pixel_percentage * 100.0,
                VALID_PIX_THRESHOLD * 100.0
            ));
        }
        if q.seamline_quality < SEAM_THRESHOLD {
            warnings.push(format!(
                "Overlap {}-{} seamline quality {:.2} below {:.2}",
                region.swath1_id, region.swath2_id, q.seamline_quality, SEAM_THRESHOLD
            ));
        }
    }

    let overall_quality = min_phase.min(min_radiometric).min(min_valid).min(min_seam);
    let validation_passed = warnings.is_empty() && overall_quality >= PHASE_THRESHOLD;

    Ok(QualityResults {
        overall_quality,
        phase_preservation: Some(min_phase),
        radiometric_consistency: min_radiometric,
        overlap_qualities: overlaps.iter().map(|o| o.quality_metrics.clone()).collect(),
        validation_passed,
        warnings,
    })
}
