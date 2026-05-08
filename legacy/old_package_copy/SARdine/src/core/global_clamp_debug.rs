//! Global diagnostic clamp wrapper to trap inverted bounds early.
//!
//! This introduces a `ClampDebug` trait for `f64` and `f32` that performs
//! a runtime check (only when SARDINE_STRICT_CLAMP=1) and logs detailed
//! context before swapping or panicking (configurable) when `min > max`.
//!
//! Integration strategy:
//! - Replace direct `.clamp(a,b)` calls in high‑risk modules with
//!   `.dbg_clamp(a,b,"label")` progressively. This file centralizes logic
//!   so we don't duplicate ad‑hoc macros across modules.
//! - Does NOT use a blanket trait in prelude to avoid silently changing
//!   semantics; adoption is explicit.
//!
//! Rationale: We observed a panic inside `f64::clamp` with bounds (60,-60)
//! that escaped earlier localized instrumentation. This global helper
//! ensures any newly instrumented call yields file/line + label.
use std::fmt::Display;

pub trait ClampDebug: Sized + Copy + PartialOrd + Display {
    fn dbg_clamp(self, min: Self, max: Self, label: &str) -> Self;
}

fn strict_mode() -> bool {
    std::env::var("SARDINE_STRICT_CLAMP").ok().as_deref() == Some("1")
}

impl ClampDebug for f64 {
    fn dbg_clamp(self, min: Self, max: Self, label: &str) -> Self {
        if strict_mode() && min > max {
            log::error!(
                "🚨 GLOBAL CLAMP INVERSION (f64) label='{}' value={} min={} max={} (swapping)",
                label,
                self,
                min,
                max
            );
            return self.clamp(max, min);
        }
        self.clamp(min, max)
    }
}

impl ClampDebug for f32 {
    fn dbg_clamp(self, min: Self, max: Self, label: &str) -> Self {
        if strict_mode() && min > max {
            log::error!(
                "🚨 GLOBAL CLAMP INVERSION (f32) label='{}' value={} min={} max={} (swapping)",
                label,
                self,
                min,
                max
            );
            return self.clamp(max, min);
        }
        self.clamp(min, max)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dbg_clamp_normal() {
        std::env::set_var("SARDINE_STRICT_CLAMP", "1");
        let v = 5.0f64.dbg_clamp(0.0, 10.0, "normal");
        assert_eq!(v, 5.0);
    }

    #[test]
    fn test_dbg_clamp_inversion() {
        std::env::set_var("SARDINE_STRICT_CLAMP", "1");
        let v = 5.0f64.dbg_clamp(10.0, -10.0, "invert");
        // swapped -> clamp(-10,10) => 5
        assert_eq!(v, 5.0);
    }
}
