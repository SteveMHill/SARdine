"""
Comprehensive RTC Metadata Documentation Module
==============================================

This module provides comprehensive metadata documentation for Range-Doppler 
Terrain Correction (RTC) processing, ensuring scientific traceability and 
reproducibility of SAR processing results.

Scientific Standards Compliance:
- ESA Sentinel-1 Product Specification: S1-RS-MDA-52-7441
- CEOS Analysis Ready Data (ARD) requirements
- IEEE SAR processing standards
- Open Geospatial Consortium (OGC) metadata standards
"""

import json
import time
from datetime import datetime, timezone
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, asdict
from pathlib import Path

@dataclass
class DEMMetadata:
    """DEM source and processing metadata"""
    source: str                    # e.g., "SRTM_30m", "Copernicus_30m", "ASTER_30m"
    version: str                   # DEM version/date
    resolution_meters: float       # Native DEM resolution in meters
    vertical_datum: str           # e.g., "EGM96", "WGS84"
    horizontal_datum: str         # e.g., "WGS84"
    coverage_area: Tuple[float, float, float, float]  # [min_lon, min_lat, max_lon, max_lat]
    void_fill_method: Optional[str] = None  # e.g., "interpolation", "nearest_neighbor"
    smoothing_applied: bool = False
    quality_flags: Dict[str, Any] = None

@dataclass
class RTCProcessingMetadata:
    """RTC method and processing parameters metadata"""
    rtc_method: str               # e.g., "gamma_flattening", "sigma_flattening", "area_projection"
    correction_formula: str       # Mathematical formula used (e.g., "γ⁰ = σ⁰ / cos(θ_local)")
    local_incidence_method: str   # e.g., "DEM_normals", "orbital_geometry", "annotation_based"
    cosine_clip_threshold: float  # ε threshold for cos(θ_local) clipping
    masking_applied: bool         # Whether quality masking was applied
    shadow_layover_detection: bool # Whether shadow/layover detection was performed
    reference_incidence_angle: Optional[float] = None  # Reference angle for some methods
    processing_mode: str = "linear"  # "linear" or "dB" processing
    
@dataclass
class CalibrationMetadata:
    """Radiometric calibration metadata"""
    calibration_type: str         # e.g., "sigma0", "beta0", "gamma0"
    calibration_equation: str     # ESA equation reference
    lut_source: str              # Path to calibration LUT XML
    noise_removal_applied: bool   # Whether thermal noise removal was applied
    noise_source: str            # Path to noise XML file
    antenna_pattern_correction: bool = False
    range_spreading_correction: bool = False

@dataclass
class ProcessingChainMetadata:
    """Complete processing chain documentation"""
    processing_steps: List[Dict[str, Any]]  # Ordered list of processing steps
    software_version: str         # SARdine version
    processing_timestamp: str     # ISO format timestamp
    input_product: str           # Original Sentinel-1 product ID
    orbit_files_applied: List[str] # List of precise orbit files used
    total_processing_time: float  # Total processing time in seconds
    
@dataclass
class QualityMetadata:
    """Quality assessment and validation metadata"""
    valid_pixel_percentage: float
    shadow_pixel_percentage: float
    layover_pixel_percentage: float
    steep_terrain_percentage: float  # Pixels with cos(θ) < threshold
    backscatter_statistics: Dict[str, float]  # min, max, mean, std
    incidence_angle_range: Tuple[float, float]  # [min, max] in degrees
    quality_flags: List[str]      # Quality warnings/flags
    
@dataclass
class GeospatialMetadata:
    """Geospatial reference information"""
    coordinate_system: str        # e.g., "EPSG:4326", "EPSG:32633"
    pixel_spacing: Tuple[float, float]  # [range, azimuth] in meters
    scene_center_lat: float
    scene_center_lon: float
    scene_extent: Tuple[float, float, float, float]  # [min_lon, min_lat, max_lon, max_lat]
    acquisition_geometry: Dict[str, Any]  # Look direction, incidence angles, etc.

@dataclass
class ComprehensiveRTCMetadata:
    """Complete RTC processing metadata package"""
    product_info: Dict[str, str]           # Product ID, polarization, etc.
    dem_metadata: DEMMetadata
    rtc_processing: RTCProcessingMetadata
    calibration_metadata: CalibrationMetadata
    processing_chain: ProcessingChainMetadata
    quality_metadata: QualityMetadata
    geospatial_metadata: GeospatialMetadata
    
    # Standard metadata fields
    creation_date: str
    format_version: str = "1.0"
    metadata_standard: str = "SARdine_RTC_v1.0"
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization"""
        return asdict(self)
    
    def to_json(self, indent: int = 2) -> str:
        """Convert to JSON string"""
        return json.dumps(self.to_dict(), indent=indent, default=str)
    
    def save_to_file(self, filepath: str) -> None:
        """Save metadata to JSON file"""
        with open(filepath, 'w') as f:
            f.write(self.to_json())
    
    @classmethod
    def from_json(cls, json_str: str) -> 'ComprehensiveRTCMetadata':
        """Create from JSON string"""
        data = json.loads(json_str)
        
        # Reconstruct nested dataclasses
        data['dem_metadata'] = DEMMetadata(**data['dem_metadata'])
        data['rtc_processing'] = RTCProcessingMetadata(**data['rtc_processing'])
        data['calibration_metadata'] = CalibrationMetadata(**data['calibration_metadata'])
        data['processing_chain'] = ProcessingChainMetadata(**data['processing_chain'])
        data['quality_metadata'] = QualityMetadata(**data['quality_metadata'])
        data['geospatial_metadata'] = GeospatialMetadata(**data['geospatial_metadata'])
        
        return cls(**data)
    
    @classmethod
    def from_file(cls, filepath: str) -> 'ComprehensiveRTCMetadata':
        """Load metadata from JSON file"""
        with open(filepath, 'r') as f:
            return cls.from_json(f.read())


class RTCMetadataBuilder:
    """Builder class for creating comprehensive RTC metadata"""
    
    def __init__(self, product_id: str, polarization: str):
        self.product_id = product_id
        self.polarization = polarization
        self.processing_steps = []
        self.start_time = time.time()
        
    def add_processing_step(self, step_name: str, step_number: int, 
                          parameters: Dict[str, Any], duration: float,
                          input_shape: Optional[Tuple[int, int]] = None,
                          output_shape: Optional[Tuple[int, int]] = None) -> None:
        """Add a processing step to the metadata"""
        step = {
            'step_number': step_number,
            'step_name': step_name,
            'parameters': parameters,
            'processing_time_seconds': duration,
            'timestamp': datetime.now(timezone.utc).isoformat(),
            'input_shape': input_shape,
            'output_shape': output_shape
        }
        self.processing_steps.append(step)
    
    def build_metadata(self, 
                      # DEM parameters
                      dem_source: str, dem_resolution: float, dem_version: str,
                      dem_coverage: Tuple[float, float, float, float],
                      # RTC parameters  
                      rtc_method: str, cosine_threshold: float, masking_enabled: bool,
                      # Calibration parameters
                      calibration_type: str, lut_source: str, noise_applied: bool,
                      # Quality metrics
                      valid_pixels_pct: float, backscatter_stats: Dict[str, float],
                      incidence_range: Tuple[float, float],
                      # Geospatial info
                      pixel_spacing: Tuple[float, float], scene_center: Tuple[float, float],
                      scene_extent: Tuple[float, float, float, float],
                      coordinate_system: str = "EPSG:4326") -> ComprehensiveRTCMetadata:
        """Build complete RTC metadata"""
        
        total_time = time.time() - self.start_time
        
        # Build DEM metadata
        dem_metadata = DEMMetadata(
            source=dem_source,
            version=dem_version,
            resolution_meters=dem_resolution,
            vertical_datum="EGM96",  # Standard for SRTM/Copernicus
            horizontal_datum="WGS84",
            coverage_area=dem_coverage,
            void_fill_method="interpolation" if any("fill" in step['step_name'] for step in self.processing_steps) else None,
            smoothing_applied=any("smooth" in step['step_name'] for step in self.processing_steps),
            quality_flags={}
        )
        
        # Build RTC processing metadata
        rtc_processing = RTCProcessingMetadata(
            rtc_method=rtc_method,
            correction_formula="γ⁰ = σ⁰ / cos(θ_local)" if "gamma" in rtc_method else "custom",
            local_incidence_method="DEM_normals",
            cosine_clip_threshold=cosine_threshold,
            masking_applied=masking_enabled,
            shadow_layover_detection=masking_enabled,
            processing_mode="linear"
        )
        
        # Build calibration metadata
        calibration_metadata = CalibrationMetadata(
            calibration_type=calibration_type,
            calibration_equation="ESA_S1_TN_MDA_52_7448_v2.0_corrected",
            lut_source=lut_source,
            noise_removal_applied=noise_applied,
            noise_source="annotation_derived"
        )
        
        # Build processing chain metadata
        processing_chain = ProcessingChainMetadata(
            processing_steps=self.processing_steps,
            software_version="SARdine v1.0",  # Should be extracted from package
            processing_timestamp=datetime.now(timezone.utc).isoformat(),
            input_product=self.product_id,
            orbit_files_applied=["precise_orbit_extracted"],  # Should be actual orbit files
            total_processing_time=total_time
        )
        
        # Build quality metadata
        quality_metadata = QualityMetadata(
            valid_pixel_percentage=valid_pixels_pct,
            shadow_pixel_percentage=100.0 - valid_pixels_pct,  # Approximation
            layover_pixel_percentage=0.0,  # Should be calculated
            steep_terrain_percentage=0.0,  # Should be calculated
            backscatter_statistics=backscatter_stats,
            incidence_angle_range=incidence_range,
            quality_flags=[]
        )
        
        # Build geospatial metadata
        geospatial_metadata = GeospatialMetadata(
            coordinate_system=coordinate_system,
            pixel_spacing=pixel_spacing,
            scene_center_lat=scene_center[1],
            scene_center_lon=scene_center[0],
            scene_extent=scene_extent,
            acquisition_geometry={
                'look_direction': 'right',
                'polarization': self.polarization,
                'incidence_angle_range': incidence_range
            }
        )
        
        # Create comprehensive metadata
        return ComprehensiveRTCMetadata(
            product_info={
                'product_id': self.product_id,
                'polarization': self.polarization,
                'processing_level': 'RTC',
                'processor': 'SARdine'
            },
            dem_metadata=dem_metadata,
            rtc_processing=rtc_processing,
            calibration_metadata=calibration_metadata,
            processing_chain=processing_chain,
            quality_metadata=quality_metadata,
            geospatial_metadata=geospatial_metadata,
            creation_date=datetime.now(timezone.utc).isoformat(),
            format_version="1.0",
            metadata_standard="SARdine_RTC_v1.0"
        )


def create_rtc_metadata_template() -> Dict[str, Any]:
    """Create a template for RTC metadata with documentation"""
    return {
        "metadata_standard": "SARdine_RTC_v1.0",
        "description": "Comprehensive metadata for SAR Range-Doppler Terrain Correction processing",
        "required_fields": {
            "dem_metadata": {
                "source": "DEM data source (e.g., SRTM_30m, Copernicus_30m)",
                "resolution_meters": "Native DEM resolution in meters",
                "vertical_datum": "Vertical reference system (e.g., EGM96)",
                "horizontal_datum": "Horizontal reference system (e.g., WGS84)"
            },
            "rtc_processing": {
                "rtc_method": "RTC algorithm used (e.g., gamma_flattening)",
                "correction_formula": "Mathematical formula applied",
                "cosine_clip_threshold": "Threshold for cosine clipping (ε)",
                "masking_applied": "Whether quality masking was performed"
            },
            "calibration_metadata": {
                "calibration_type": "Radiometric calibration type (sigma0, beta0, gamma0)",
                "calibration_equation": "ESA equation reference",
                "noise_removal_applied": "Whether thermal noise removal was applied"
            },
            "quality_metadata": {
                "valid_pixel_percentage": "Percentage of valid pixels",
                "backscatter_statistics": "Statistical summary of backscatter values",
                "incidence_angle_range": "Range of local incidence angles"
            }
        },
        "scientific_references": [
            "ESA Sentinel-1 Product Specification: S1-RS-MDA-52-7441",
            "Small, D. (2011): Flattening gamma: Radiometric terrain correction for SAR imagery",
            "Ulaby & Long (2014): Microwave Radar and Radiometric Remote Sensing"
        ]
    }


# Example usage and testing
if __name__ == "__main__":
    # Create example metadata
    builder = RTCMetadataBuilder("S1A_IW_20240720T052639_DVP_RTC20_G", "VV")
    
    # Add some example processing steps
    builder.add_processing_step("Metadata Extraction", 1, {"source": "annotation_xml"}, 0.5)
    builder.add_processing_step("Precise Orbit Application", 2, {"orbit_type": "precise"}, 2.1)
    builder.add_processing_step("IW Merging", 3, {"subswaths": ["IW1", "IW2", "IW3"]}, 5.2)
    builder.add_processing_step("Terrain Correction", 4, {"method": "gamma_flattening", "cosine_threshold": 0.1}, 15.8)
    
    # Build comprehensive metadata
    metadata = builder.build_metadata(
        # DEM parameters
        dem_source="SRTM_30m",
        dem_resolution=30.0,
        dem_version="v3.0",
        dem_coverage=(-122.5, 37.5, -122.0, 38.0),
        # RTC parameters
        rtc_method="gamma_flattening_with_masking",
        cosine_threshold=0.1,
        masking_enabled=True,
        # Calibration parameters
        calibration_type="sigma0",
        lut_source="calibration.xml",
        noise_applied=True,
        # Quality metrics
        valid_pixels_pct=95.2,
        backscatter_stats={"min": -30.5, "max": 15.2, "mean": -12.8, "std": 4.3},
        incidence_range=(29.1, 45.9),
        # Geospatial info
        pixel_spacing=(20.0, 20.0),
        scene_center=(-122.25, 37.75),
        scene_extent=(-122.5, 37.5, -122.0, 38.0)
    )
    
    # Print metadata
    print("🚀 Comprehensive RTC Metadata Example")
    print("=" * 50)
    print(metadata.to_json())
    
    # Save to file
    output_file = "/tmp/rtc_metadata_example.json"
    metadata.save_to_file(output_file)
    print(f"\n✅ Metadata saved to: {output_file}")
    
    # Test loading
    loaded_metadata = ComprehensiveRTCMetadata.from_file(output_file)
    print(f"\n✅ Successfully loaded metadata with {len(loaded_metadata.processing_chain.processing_steps)} processing steps")