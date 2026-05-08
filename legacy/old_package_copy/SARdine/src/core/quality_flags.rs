//! Quality flags for scientific processing validation.
//!
//! This module implements the quality flag system required by the scientific audit
//! (NOISE-1, DEBURST-1, RTC-1). Quality flags are surfaced in output metadata
//! instead of silent fallbacks.
//!
//! # Usage
//!
//! ```rust,ignore
//! use sardine::core::quality_flags::{QualityFlags, QualityFlag};
//!
//! let mut flags = QualityFlags::new();
//! flags.add(QualityFlag::NoiseRemovalSkipped {
//!     reason: "Zero ratio exceeded threshold".to_string(),
//!     zero_ratio: 0.96,
//!     threshold: 0.95,
//! });
//!
//! assert!(flags.has_flag("NOISE_REMOVAL_SKIPPED"));
//! ```

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::sync::{Arc, RwLock};

/// Individual quality flag with detailed context.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "flag_type")]
pub enum QualityFlag {
    /// NOISE-1: Thermal noise removal was skipped due to threshold checks.
    NoiseRemovalSkipped {
        reason: String,
        zero_ratio: Option<f64>,
        noise_power_ratio: Option<f64>,
        zero_threshold: f64,
        ratio_threshold: f64,
    },

    /// NOISE-1: Noise removal thresholds were modified by environment variables.
    NoiseThresholdsModified {
        zero_threshold: f64,
        ratio_threshold: f64,
        default_zero: f64,
        default_ratio: f64,
    },

    /// DEBURST-1: Polynomial t0 fallback occurred.
    DeburstT0Fallback {
        burst_id: usize,
        reason: String,
        expected_offset: Option<f64>,
    },

    /// DEBURST-1: Time domain mismatch detected in deburst.
    DeburstTimeDomainMismatch {
        burst_id: usize,
        burst_time: f64,
        poly_t0: f64,
        offset: f64,
    },

    /// RTC-1: Default incidence angle was used instead of annotation value.
    RtcDefaultIncidenceAngle {
        default_angle_deg: f64,
        reason: String,
    },

    /// RTC: Incidence angle profile was missing or incomplete.
    RtcIncidenceProfileMissing { reason: String },

    /// Calibration LUT parsing issue.
    CalibrationLutWarning { lut_type: String, issue: String },

    /// Generic warning flag for other conditions.
    Warning {
        code: String,
        message: String,
        details: HashMap<String, String>,
    },
}

impl QualityFlag {
    /// Get the flag type name for lookup.
    pub fn flag_type(&self) -> &'static str {
        match self {
            QualityFlag::NoiseRemovalSkipped { .. } => "NOISE_REMOVAL_SKIPPED",
            QualityFlag::NoiseThresholdsModified { .. } => "NOISE_THRESHOLDS_MODIFIED",
            QualityFlag::DeburstT0Fallback { .. } => "DEBURST_T0_FALLBACK",
            QualityFlag::DeburstTimeDomainMismatch { .. } => "DEBURST_TIME_DOMAIN_MISMATCH",
            QualityFlag::RtcDefaultIncidenceAngle { .. } => "RTC_DEFAULT_INCIDENCE_ANGLE",
            QualityFlag::RtcIncidenceProfileMissing { .. } => "RTC_INCIDENCE_PROFILE_MISSING",
            QualityFlag::CalibrationLutWarning { .. } => "CALIBRATION_LUT_WARNING",
            QualityFlag::Warning {  .. } => {
                // This is a bit awkward but we need a static str
                // In practice, use the code field for lookup
                "WARNING"
            }
        }
    }

    /// Get the severity level (for filtering/reporting).
    pub fn severity(&self) -> &'static str {
        match self {
            QualityFlag::NoiseRemovalSkipped { .. } => "WARNING",
            QualityFlag::NoiseThresholdsModified { .. } => "INFO",
            QualityFlag::DeburstT0Fallback { .. } => "WARNING",
            QualityFlag::DeburstTimeDomainMismatch { .. } => "ERROR",
            QualityFlag::RtcDefaultIncidenceAngle { .. } => "WARNING",
            QualityFlag::RtcIncidenceProfileMissing { .. } => "WARNING",
            QualityFlag::CalibrationLutWarning { .. } => "WARNING",
            QualityFlag::Warning { .. } => "WARNING",
        }
    }
}

/// Collection of quality flags for a processing run.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct QualityFlags {
    flags: Vec<QualityFlag>,
}

impl QualityFlags {
    /// Create an empty quality flags collection.
    pub fn new() -> Self {
        Self { flags: Vec::new() }
    }

    /// Add a quality flag.
    pub fn add(&mut self, flag: QualityFlag) {
        log::info!("🚩 Quality flag: {} - {:?}", flag.flag_type(), flag);
        self.flags.push(flag);
    }

    /// Check if a specific flag type is present.
    pub fn has_flag(&self, flag_type: &str) -> bool {
        self.flags.iter().any(|f| f.flag_type() == flag_type)
    }

    /// Get all flags.
    pub fn all_flags(&self) -> &[QualityFlag] {
        &self.flags
    }

    /// Get flags by severity.
    pub fn flags_by_severity(&self, severity: &str) -> Vec<&QualityFlag> {
        self.flags
            .iter()
            .filter(|f| f.severity() == severity)
            .collect()
    }

    /// Check if any error-level flags are present.
    pub fn has_errors(&self) -> bool {
        self.flags.iter().any(|f| f.severity() == "ERROR")
    }

    /// Check if any warning-level flags are present.
    pub fn has_warnings(&self) -> bool {
        self.flags.iter().any(|f| f.severity() == "WARNING")
    }

    /// Get count of flags.
    pub fn count(&self) -> usize {
        self.flags.len()
    }

    /// Merge flags from another collection.
    pub fn merge(&mut self, other: QualityFlags) {
        self.flags.extend(other.flags);
    }

    /// Convert to JSON string for output metadata.
    pub fn to_json(&self) -> String {
        serde_json::to_string_pretty(&self).unwrap_or_else(|_| "[]".to_string())
    }

    /// Convert to a summary map for embedding in metadata.
    pub fn to_summary(&self) -> HashMap<String, Vec<String>> {
        let mut summary: HashMap<String, Vec<String>> = HashMap::new();

        for flag in &self.flags {
            let key = flag.flag_type().to_string();
            let msg = match flag {
                QualityFlag::NoiseRemovalSkipped { reason, .. } => reason.clone(),
                QualityFlag::NoiseThresholdsModified {
                    zero_threshold,
                    ratio_threshold,
                    ..
                } => {
                    format!("zero={:.2}, ratio={:.2}", zero_threshold, ratio_threshold)
                }
                QualityFlag::DeburstT0Fallback {
                    burst_id, reason, ..
                } => {
                    format!("burst {}: {}", burst_id, reason)
                }
                QualityFlag::DeburstTimeDomainMismatch {
                    burst_id, offset, ..
                } => {
                    format!("burst {}: offset={:.1}s", burst_id, offset)
                }
                QualityFlag::RtcDefaultIncidenceAngle {
                    default_angle_deg,
                    reason,
                } => {
                    format!("{:.1}° used: {}", default_angle_deg, reason)
                }
                QualityFlag::RtcIncidenceProfileMissing { reason } => reason.clone(),
                QualityFlag::CalibrationLutWarning { lut_type, issue } => {
                    format!("{}: {}", lut_type, issue)
                }
                QualityFlag::Warning { code, message, .. } => {
                    format!("{}: {}", code, message)
                }
            };
            summary.entry(key).or_insert_with(Vec::new).push(msg);
        }

        summary
    }
}

/// Thread-safe global quality flag collector.
/// Used to accumulate flags across parallel processing stages.
#[derive(Clone)]
pub struct GlobalQualityFlags {
    inner: Arc<RwLock<QualityFlags>>,
}

impl GlobalQualityFlags {
    /// Create a new global quality flag collector.
    pub fn new() -> Self {
        Self {
            inner: Arc::new(RwLock::new(QualityFlags::new())),
        }
    }

    /// Add a flag (thread-safe).
    pub fn add(&self, flag: QualityFlag) {
        if let Ok(mut flags) = self.inner.write() {
            flags.add(flag);
        }
    }

    /// Get a copy of all flags.
    pub fn get_flags(&self) -> QualityFlags {
        self.inner.read().map(|f| f.clone()).unwrap_or_default()
    }

    /// Clear all flags.
    pub fn clear(&self) {
        if let Ok(mut flags) = self.inner.write() {
            flags.flags.clear();
        }
    }
}

impl Default for GlobalQualityFlags {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// GLOBAL QUALITY FLAG SINGLETON
// ============================================================================

use std::sync::OnceLock;

static GLOBAL_QUALITY_FLAGS: OnceLock<GlobalQualityFlags> = OnceLock::new();

/// Get the global quality flag collector (singleton).
///
/// This allows processing functions (noise removal, deburst, RTC) to emit
/// quality flags without requiring API changes. The flags can be retrieved
/// at the end of processing for inclusion in output metadata.
///
/// # Example
///
/// ```rust,ignore
/// use sardine::core::quality_flags::{global_quality_flags, QualityFlag};
///
/// // During processing:
/// global_quality_flags().add(QualityFlag::NoiseRemovalSkipped { ... });
///
/// // At end of processing:
/// let flags = global_quality_flags().get_flags();
/// let json = flags.to_json();
/// ```
pub fn global_quality_flags() -> &'static GlobalQualityFlags {
    GLOBAL_QUALITY_FLAGS.get_or_init(GlobalQualityFlags::new)
}

/// Reset the global quality flags (for testing or between processing runs).
/// Note: This creates a new GlobalQualityFlags, the old one becomes inaccessible.
pub fn reset_global_quality_flags() {
    if let Some(flags) = GLOBAL_QUALITY_FLAGS.get() {
        flags.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quality_flag_creation() {
        let flag = QualityFlag::NoiseRemovalSkipped {
            reason: "Test reason".to_string(),
            zero_ratio: Some(0.96),
            noise_power_ratio: Some(0.8),
            zero_threshold: 0.95,
            ratio_threshold: 0.9,
        };

        assert_eq!(flag.flag_type(), "NOISE_REMOVAL_SKIPPED");
        assert_eq!(flag.severity(), "WARNING");
    }

    #[test]
    fn test_quality_flags_collection() {
        let mut flags = QualityFlags::new();

        flags.add(QualityFlag::DeburstT0Fallback {
            burst_id: 1,
            reason: "Missing t0".to_string(),
            expected_offset: None,
        });

        assert!(flags.has_flag("DEBURST_T0_FALLBACK"));
        assert!(!flags.has_flag("NOISE_REMOVAL_SKIPPED"));
        assert_eq!(flags.count(), 1);
    }

    #[test]
    fn test_quality_flags_json() {
        let mut flags = QualityFlags::new();
        flags.add(QualityFlag::RtcDefaultIncidenceAngle {
            default_angle_deg: 35.0,
            reason: "Annotation not available".to_string(),
        });

        let json = flags.to_json();
        // The flag_type is serialized as a tag by serde
        assert!(json.contains("RtcDefaultIncidenceAngle") || json.contains("flag_type"));
        assert!(json.contains("35.0") || json.contains("35"));
    }

    #[test]
    fn test_global_quality_flags_thread_safe() {
        let global = GlobalQualityFlags::new();

        let g1 = global.clone();
        let g2 = global.clone();

        std::thread::spawn(move || {
            g1.add(QualityFlag::Warning {
                code: "TEST1".to_string(),
                message: "Thread 1".to_string(),
                details: HashMap::new(),
            });
        })
        .join()
        .unwrap();

        std::thread::spawn(move || {
            g2.add(QualityFlag::Warning {
                code: "TEST2".to_string(),
                message: "Thread 2".to_string(),
                details: HashMap::new(),
            });
        })
        .join()
        .unwrap();

        let flags = global.get_flags();
        assert_eq!(flags.count(), 2);
    }
}
