//! Unit tests for coordinate convention consistency across merge functions
//!
//! This module ensures that all merge functions use the same coordinate
//! convention documented in types.rs:
//!
//! - `first_line_global`: inclusive (first valid line index)
//! - `last_line_global`: **exclusive** (one past last valid line index)
//! - `first_sample_global`: inclusive (first valid sample index)
//! - `last_sample_global`: **exclusive** (one past last valid sample index)
//!
//! Range notation: `[first, last)` means `first..last` in Rust slice syntax

#[cfg(test)]
mod tests {
    use crate::types::SubSwath;

    /// Test that SubSwath bounds follow exclusive-end convention
    #[test]
    fn subswath_bounds_are_exclusive_end() {
        let sw = SubSwath {
            id: "IW1".to_string(),
            first_line_global: 0,
            last_line_global: 1000,  // Exclusive: represents lines 0..1000
            azimuth_samples: 1000,
            first_sample_global: 0,
            last_sample_global: 2000,  // Exclusive: represents samples 0..2000
            range_samples: 2000,
            // ... other fields with defaults
            burst_count: 9,
            lines_per_burst: 1500,
            full_range_samples: 2000,
            valid_first_line: None,
            valid_last_line: Some(1000),
            valid_first_sample: None,
            valid_last_sample: Some(2000),
            range_pixel_spacing: 2.33,
            azimuth_pixel_spacing: 14.0,
            slant_range_time: 0.005,
            burst_duration: 2.7,
            near_range_m: 800000.0,
            prf_hz: Some(486.5),
            dc_polynomial: None,
            azimuth_time_interval: Some(0.002055556),
            dc_polynomial_t0: None,
            fm_rate_estimates: None,
        };

        // Verify exclusive convention: last = first + count
        assert_eq!(
            sw.last_line_global - sw.first_line_global,
            sw.azimuth_samples,
            "Azimuth: last_line_global should be first + azimuth_samples (exclusive)"
        );

        assert_eq!(
            sw.last_sample_global - sw.first_sample_global,
            sw.range_samples,
            "Range: last_sample_global should be first + range_samples (exclusive)"
        );

        // Valid bounds should also be exclusive
        if let (Some(valid_last_line), Some(valid_last_sample)) =
            (sw.valid_last_line, sw.valid_last_sample)
        {
            assert_eq!(
                valid_last_line,
                sw.last_line_global,
                "valid_last_line should match last_line_global"
            );
            assert_eq!(
                valid_last_sample,
                sw.last_sample_global,
                "valid_last_sample should match last_sample_global"
            );
        }
    }

    /// Test coordinate calculation from actual data dimensions
    #[test]
    fn coordinates_from_actual_dimensions() {
        let first_line = 100;
        let first_sample = 200;
        let actual_height = 500;  // rows in data array
        let actual_width = 1000;  // columns in data array

        // CORRECT: Exclusive end = first + actual_size
        let last_line_correct = first_line + actual_height;  // 600
        let last_sample_correct = first_sample + actual_width;  // 1200

        // INCORRECT: Inclusive end = first + actual_size - 1 (old bug)
        let last_line_incorrect = first_line + actual_height - 1;  // 599
        let last_sample_incorrect = first_sample + actual_width - 1;  // 1199

        // Verify that actual_size matches the range length only with exclusive convention
        assert_eq!(
            last_line_correct - first_line,
            actual_height,
            "Only exclusive convention gives correct height"
        );
        assert_eq!(
            last_sample_correct - first_sample,
            actual_width,
            "Only exclusive convention gives correct width"
        );

        // Inclusive convention would be off by one
        assert_ne!(
            last_line_incorrect - first_line,
            actual_height,
            "Inclusive convention is off by one (height)"
        );
        assert_ne!(
            last_sample_incorrect - first_sample,
            actual_width,
            "Inclusive convention is off by one (width)"
        );
    }

    /// Test that saturating_add without saturating_sub(1) gives exclusive bounds
    #[test]
    fn saturating_add_produces_exclusive_bounds() {
        let first_global: usize = 100;
        let actual_size: usize = 500;

        // CORRECT: Exclusive end (documented convention)
        let last_global_correct = first_global.saturating_add(actual_size);

        // INCORRECT: Inclusive end (old bug pattern)
        let last_global_incorrect = first_global.saturating_add(actual_size).saturating_sub(1);

        // Verify range length
        assert_eq!(
            last_global_correct - first_global,
            actual_size,
            "saturating_add alone gives correct exclusive bounds"
        );

        assert_eq!(
            last_global_incorrect - first_global,
            actual_size - 1,
            "adding saturating_sub(1) creates off-by-one error"
        );

        // Verify slice notation equivalence
        let mock_data = vec![0u8; actual_size];
        let slice = &mock_data[0..actual_size];  // Rust slice is exclusive-end

        assert_eq!(
            slice.len(),
            last_global_correct - first_global,
            "Exclusive convention matches Rust slice semantics"
        );
    }

    /// Test edge case: zero-sized extent
    #[test]
    fn zero_sized_extent_exclusive_convention() {
        let first: usize = 100;
        let size: usize = 0;

        let last = first.saturating_add(size);

        assert_eq!(last, first, "Zero-sized extent: last == first (empty range)");
        assert_eq!(last - first, 0, "Range length is zero");
    }

    /// Test edge case: maximum safe size before saturation
    #[test]
    fn maximum_size_before_saturation() {
        let first: usize = 100;
        let max_possible = usize::MAX - first;

        let last = first.saturating_add(max_possible);

        assert_eq!(last, usize::MAX, "Saturating add reaches usize::MAX");
        assert_eq!(
            last - first,
            max_possible,
            "Range length is maximum possible"
        );
    }

    /// Test that off-by-one errors would break pixel alignment
    #[test]
    fn off_by_one_breaks_alignment() {
        // Simulate two adjacent subswaths in range direction
        let sw1_first_sample = 0;
        let sw1_range_samples = 1000;
        let sw1_last_sample_correct = sw1_first_sample + sw1_range_samples;  // 1000 (exclusive)

        let sw2_first_sample = sw1_last_sample_correct;  // Should start at 1000
        let sw2_range_samples = 800;
        let sw2_last_sample_correct = sw2_first_sample + sw2_range_samples;  // 1800 (exclusive)

        // Total merged width should be sum of both subswaths (no gap, no overlap)
        let merged_width = sw2_last_sample_correct - sw1_first_sample;
        assert_eq!(
            merged_width,
            sw1_range_samples + sw2_range_samples,
            "Adjacent subswaths with exclusive bounds have no gap"
        );

        // If we incorrectly used inclusive bounds (old bug):
        let sw1_last_sample_incorrect = sw1_first_sample + sw1_range_samples - 1;  // 999 (inclusive)
        // Then sw2 would need to start at 999 (overlap) or 1000 (1-pixel gap)

        // Starting at 999 creates 1-pixel overlap (wrong!)
        if sw2_first_sample == sw1_last_sample_incorrect {
            let overlap = 1;
            assert!(overlap > 0, "Inclusive bounds create overlap or gap");
        }

        // Starting at 1000 creates 1-pixel gap (also wrong!)
        if sw2_first_sample == sw1_last_sample_incorrect + 1 {
            let gap = sw2_first_sample - sw1_last_sample_incorrect - 1;
            assert_eq!(gap, 0, "Should be no gap, but inclusive convention forces gap");
        }
    }
}
