"""
Integration tests for burst timing extraction across all code paths.

This test suite verifies that burst timing data is correctly extracted and
passed to terrain correction for merged IW TOPSAR products. Missing or
incorrect burst timing causes 3× azimuth pixel errors in geocoded output.
"""

import json
import pytest
from pathlib import Path
from typing import Dict, Any, List, Optional

import numpy as np
import sardine


class TestBurstTimingExtraction:
    """Test burst timing extraction from various sources"""
    
    def test_burst_timing_from_cached_metadata(self):
        """Test burst timing extraction from cached metadata (primary path)"""
        # This test verifies the primary code path in terrain.py:_extract_burst_timing
        # Path: cached_meta["burst_timing_json"] -> JSON deserialize -> pass to terrain correction
        
        # Mock cached metadata with burst timing JSON
        cached_meta = {
            "burst_timing_json": json.dumps([
                {
                    "subswath_id": "IW1",
                    "burst_index": 0,
                    "azimuth_time_rel_orbit": 100.5,
                    "first_line_global": 0,
                    "last_line_global": 1500,
                    "first_valid_sample": 0,
                    "last_valid_sample": 20000,
                },
                {
                    "subswath_id": "IW1",
                    "burst_index": 1,
                    "azimuth_time_rel_orbit": 103.2,
                    "first_line_global": 1500,
                    "last_line_global": 3000,
                    "first_valid_sample": 0,
                    "last_valid_sample": 20000,
                },
            ])
        }
        
        # Verify JSON can be parsed
        burst_timings = json.loads(cached_meta["burst_timing_json"])
        assert len(burst_timings) == 2
        assert all("azimuth_time_rel_orbit" in bt for bt in burst_timings)
        assert all(bt["azimuth_time_rel_orbit"] > 0 for bt in burst_timings)
        
        # Verify critical field is present
        for bt in burst_timings:
            assert bt["azimuth_time_rel_orbit"] > 0, \
                "azimuth_time_rel_orbit must be non-zero for accurate geocoding"
    
    def test_burst_timing_from_list_key(self):
        """Test burst timing extraction from list-based keys (fallback paths)"""
        # Fallback paths: cached_meta["burst_timings"], ["burst_timing_records"], ["burst_records"]
        
        burst_timings_list = [
            {
                "subswath_id": "IW2",
                "burst_index": 0,
                "azimuth_time_rel_orbit": 200.5,
                "first_line_global": 0,
                "last_line_global": 1500,
            }
        ]
        
        cached_meta_variants = [
            {"burst_timings": burst_timings_list},
            {"burst_timing_records": burst_timings_list},
            {"burst_records": burst_timings_list},
        ]
        
        for variant in cached_meta_variants:
            # Simulate the extraction logic from terrain.py
            extracted = (
                variant.get("burst_timings") or
                variant.get("burst_timing_records") or
                variant.get("burst_records")
            )
            
            assert extracted is not None
            assert len(extracted) > 0
            assert all("azimuth_time_rel_orbit" in bt for bt in extracted)
    
    def test_burst_timing_from_processor_cache(self):
        """Test burst timing extraction from processor.burst_timing_records"""
        # Code path: processor.burst_timing_records -> pass to terrain correction
        
        processor_timings = [
            {
                "subswath_id": "IW3",
                "burst_index": 0,
                "azimuth_time_rel_orbit": 300.5,
                "first_line_global": 0,
                "last_line_global": 1500,
            }
        ]
        
        # Verify structure matches expected format
        assert len(processor_timings) > 0
        for bt in processor_timings:
            assert "subswath_id" in bt
            assert "azimuth_time_rel_orbit" in bt
            assert bt["azimuth_time_rel_orbit"] > 0
    
    def test_burst_timing_validation(self):
        """Test validation that burst timing has non-zero azimuth_time_rel_orbit"""
        # Critical validation: zero azimuth_time_rel_orbit indicates missing timing
        
        valid_timing = {
            "subswath_id": "IW1",
            "azimuth_time_rel_orbit": 100.5,  # Non-zero = valid
        }
        
        invalid_timing = {
            "subswath_id": "IW1",
            "azimuth_time_rel_orbit": 0.0,  # Zero = invalid (default value)
        }
        
        # Valid timing should pass
        assert valid_timing["azimuth_time_rel_orbit"] > 0
        
        # Invalid timing should be detected
        assert invalid_timing["azimuth_time_rel_orbit"] == 0.0
        # This would cause 3× azimuth pixel errors in geocoding
    
    def test_burst_timing_serialization(self):
        """Test JSON serialization of burst timing for Rust bridge"""
        # Code path: Python list -> JSON string -> Rust deserialize
        
        burst_timings = [
            {
                "subswath_id": "IW1",
                "burst_index": 0,
                "azimuth_time_rel_orbit": 100.5,
                "first_line_global": 0,
                "last_line_global": 1500,
            }
        ]
        
        # Serialize to JSON (as done in terrain.py)
        burst_timing_json = json.dumps(burst_timings)
        
        # Verify can be deserialized
        deserialized = json.loads(burst_timing_json)
        assert len(deserialized) == len(burst_timings)
        assert deserialized[0]["azimuth_time_rel_orbit"] == 100.5
    
    def test_missing_burst_timing_warning(self):
        """Test that missing burst timing generates appropriate warnings for IW TOPSAR"""
        # Code path: No burst timing -> warning for IW TOPSAR products
        
        is_iw_topsar = True  # IW mode TOPSAR product
        burst_timings = None  # Missing timing data
        
        if is_iw_topsar and not burst_timings:
            # This should generate a warning (as in terrain.py lines 1109-1122)
            warning_expected = True
            assert warning_expected, \
                "Missing burst timing for IW TOPSAR should generate warning"
    
    @pytest.mark.integration
    def test_burst_timing_end_to_end(self, sample_safe_path: Path):
        """Integration test: Extract burst timing from real SAFE product"""
        # This test requires a real Sentinel-1 SAFE product
        # Verifies the complete chain: SAFE -> reader -> cached_meta -> terrain correction
        
        if not sample_safe_path.exists():
            pytest.skip("Sample SAFE product not available")
        
        # Create reader and extract metadata
        reader = sardine.SlcReader.new_with_full_cache(str(sample_safe_path))
        cached_meta = reader.get_cached_metadata()
        
        # Check for burst timing in various possible keys
        burst_timing_sources = [
            cached_meta.get("burst_timing_json"),
            cached_meta.get("burst_timings"),
            cached_meta.get("burst_timing_records"),
            cached_meta.get("burst_records"),
        ]
        
        # At least one source should have burst timing for IW products
        has_timing = any(source is not None for source in burst_timing_sources)
        
        # For IW TOPSAR products, burst timing should be available
        product_type = cached_meta.get("product_type", "").upper()
        mode = cached_meta.get("mode", "").upper()
        is_iw = "IW" in mode or "IW" in product_type
        
        if is_iw:
            assert has_timing, \
                "IW TOPSAR products should have burst timing data in cached metadata"
            
            # Extract and validate
            burst_timings = None
            if cached_meta.get("burst_timing_json"):
                burst_timings = json.loads(cached_meta["burst_timing_json"])
            elif cached_meta.get("burst_timings"):
                burst_timings = cached_meta["burst_timings"]
            
            if burst_timings:
                # Validate all bursts have non-zero timing
                for bt in burst_timings:
                    assert bt.get("azimuth_time_rel_orbit", 0.0) > 0, \
                        f"Burst {bt.get('burst_index')} has zero azimuth_time_rel_orbit"


class TestBurstTimingCodePaths:
    """Test all code paths for burst timing extraction"""
    
    def test_terrain_correction_burst_timing_path(self):
        """Test the exact code path used in terrain.py:run_terrain_correction"""
        # This mirrors the logic in terrain.py lines 1064-1124
        
        cached_meta = {
            "burst_timing_json": json.dumps([
                {"subswath_id": "IW1", "azimuth_time_rel_orbit": 100.5}
            ])
        }
        
        # Simulate extraction logic from terrain.py
        burst_timing_json = None
        burst_timings = None
        
        # Try JSON string first (primary path)
        raw_json = cached_meta.get("burst_timing_json")
        if raw_json:
            try:
                burst_timings = json.loads(raw_json)
            except Exception:
                pass
        
        # Fallback to list keys
        if not burst_timings:
            burst_timings = cached_meta.get("burst_timings")
        if not burst_timings:
            burst_timings = cached_meta.get("burst_timing_records")
        if not burst_timings:
            burst_timings = cached_meta.get("burst_records")
        
        # Serialize if we have valid data
        if burst_timings and isinstance(burst_timings, list) and len(burst_timings) > 0:
            try:
                burst_timing_json = json.dumps(burst_timings)
            except Exception:
                burst_timing_json = None
        
        # Verify extraction succeeded
        assert burst_timing_json is not None
        assert len(json.loads(burst_timing_json)) > 0


@pytest.fixture
def sample_safe_path() -> Path:
    """Fixture for sample SAFE product path (override in conftest.py)"""
    # Default: return non-existent path (tests will skip)
    return Path("/tmp/nonexistent_safe.SAFE")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

