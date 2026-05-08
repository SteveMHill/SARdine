"""Regression tests for the RTCMetadataBuilder construction path."""

from __future__ import annotations

from sardine.metadata.rtc_metadata import (
    RTCMetadataBuilder,
)


def test_build_metadata_produces_complete_package() -> None:
    """Ensure the builder returns a populated ComprehensiveRTCMetadata."""

    builder = RTCMetadataBuilder("S1A_TEST_PRODUCT", "VV")
    builder.add_processing_step("Test Step", 1, {"key": "value"}, duration=0.5)
    builder.set_software_version("SARdine v0.0.0-test")
    builder.add_orbit_file("S1A_TEST_ORBIT.EOF")
    builder.set_quality_metrics(
        {
            "valid_pixel_percentage": 97.5,
            "backscatter_statistics": {"min": -32.1, "max": 7.4, "mean": -12.8, "std": 4.1},
            "shadow_pixel_percentage": 1.4,
            "layover_pixel_percentage": 0.9,
            "dem_void_fill_method": "nearest_neighbor",
            "dem_smoothing_applied": True,
            "dem_quality_flags": {"void_fill": "applied"},
        }
    )

    metadata = builder.build_metadata(
        dem_source="SRTM_30m",
        dem_resolution=30.0,
        dem_version="v3",
        dem_coverage=(-122.45, 37.6, -122.0, 38.0),
        rtc_method="gamma_flattening",
        cosine_threshold=0.1,
        masking_enabled=True,
        calibration_type="sigma0",
        lut_source="calibration.xml",
        noise_applied=True,
        incidence_range=(29.0, 45.0),
        pixel_spacing=(20.0, 20.0),
        scene_center=(-122.2, 37.8),
        scene_extent=(-122.45, 37.6, -122.0, 38.0),
    )

    assert metadata.product_info["product_id"] == "S1A_TEST_PRODUCT"
    assert metadata.dem_metadata.source == "SRTM_30m"
    assert metadata.dem_metadata.coverage_area == (-122.45, 37.6, -122.0, 38.0)
    assert metadata.rtc_processing.rtc_method == "gamma_flattening"
    assert metadata.rtc_processing.reference_incidence_angle == 37.0
    assert metadata.calibration_metadata.noise_removal_applied is True
    assert metadata.quality_metadata.incidence_angle_range == (29.0, 45.0)
    assert metadata.geospatial_metadata.pixel_spacing == (20.0, 20.0)
