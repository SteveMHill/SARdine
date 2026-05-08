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
from dataclasses import dataclass, asdict, field
from pathlib import Path

# Import provenance capture - deferred to avoid circular imports
def _get_provenance() -> Dict[str, Any]:
    """Get provenance dict, importing lazily."""
    try:
        from sardine.provenance import get_provenance_dict
        return get_provenance_dict()
    except ImportError:
        return {"error": "provenance module not available"}

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
    # PROV-1/PROV-2: Full provenance including git hash and env vars
    provenance: Dict[str, Any] = field(default_factory=dict)
    
@dataclass
class QualityMetadata:
    """Quality assessment and validation metadata"""
    valid_pixel_percentage: float
    backscatter_statistics: Dict[str, float]  # min, max, mean, std
    incidence_angle_range: Tuple[float, float]  # [min, max] in degrees
    quality_flags: List[str]      # Quality warnings/flags
    shadow_pixel_percentage: Optional[float] = None
    layover_pixel_percentage: Optional[float] = None
    steep_terrain_percentage: Optional[float] = None  # Pixels with cos(θ) < threshold
    
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
        self.processing_steps: List[Dict[str, Any]] = []
        self.start_time = time.time()
        self.software_version: Optional[str] = None
        self.orbit_files: List[str] = []
        self.quality_metrics: Dict[str, Any] = {}
        
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

    def set_software_version(self, version: str) -> None:
        if not version:
            raise ValueError("Software version must be a non-empty string")
        self.software_version = version

    def add_orbit_file(self, orbit_filename: str) -> None:
        if not orbit_filename:
            raise ValueError("Orbit filename must be provided")
        self.orbit_files.append(orbit_filename)

    def set_quality_metrics(self, metrics: Dict[str, Any]) -> None:
        if not metrics:
            raise ValueError("Quality metrics dictionary cannot be empty")
        self.quality_metrics = metrics
    
    def build_metadata(self,
                       dem_source: str,
                       dem_resolution: float,
                       dem_version: str,
                       dem_coverage: Tuple[float, float, float, float],
                       rtc_method: str,
                       cosine_threshold: float,
                       masking_enabled: bool,
                       calibration_type: str,
                       lut_source: str,
                       noise_applied: bool,
                       incidence_range: Tuple[float, float],
                       pixel_spacing: Tuple[float, float],
                       scene_center: Tuple[float, float],
                       scene_extent: Tuple[float, float, float, float],
                       coordinate_system: str = "EPSG:4326") -> ComprehensiveRTCMetadata:
        """Build complete RTC metadata from verified processing artefacts."""
        
        total_time = time.time() - self.start_time

        if not self.processing_steps:
            raise ValueError("No processing steps recorded for metadata")

        if self.software_version is None:
            raise ValueError("Software version must be set via set_software_version()")

        if not self.orbit_files:
            raise ValueError("At least one precise orbit filename must be recorded")

        if 'valid_pixel_percentage' not in self.quality_metrics:
            raise ValueError("Quality metrics must include 'valid_pixel_percentage'")
        if 'backscatter_statistics' not in self.quality_metrics:
            raise ValueError("Quality metrics must include 'backscatter_statistics'")

        backscatter_stats = self.quality_metrics['backscatter_statistics']
        if not isinstance(backscatter_stats, dict) or not backscatter_stats:
            raise ValueError("backscatter_statistics must be a populated dictionary")

        shadow_pixels_pct = self.quality_metrics.get('shadow_pixel_percentage')
        layover_pixels_pct = self.quality_metrics.get('layover_pixel_percentage')
        steep_terrain_pct = self.quality_metrics.get('steep_terrain_percentage')

        try:
            dem_coverage_tuple = tuple(float(v) for v in dem_coverage)
        except (TypeError, ValueError):
            raise ValueError("dem_coverage must contain four numeric values")
        if len(dem_coverage_tuple) != 4:
            raise ValueError("dem_coverage must be a 4-element bounding box")

        try:
            pixel_spacing_tuple = tuple(float(v) for v in pixel_spacing)
        except (TypeError, ValueError):
            raise ValueError("pixel_spacing must contain two numeric values")
        if len(pixel_spacing_tuple) != 2:
            raise ValueError("pixel_spacing must be a 2-element tuple (range, azimuth)")

        try:
            scene_center_tuple = tuple(float(v) for v in scene_center)
        except (TypeError, ValueError):
            raise ValueError("scene_center must contain two numeric values")
        if len(scene_center_tuple) != 2:
            raise ValueError("scene_center must be a 2-element tuple (lon, lat)")

        try:
            scene_extent_tuple = tuple(float(v) for v in scene_extent)
        except (TypeError, ValueError):
            raise ValueError("scene_extent must contain four numeric values")
        if len(scene_extent_tuple) != 4:
            raise ValueError("scene_extent must be a 4-element bounding box")

        try:
            incidence_range_tuple = tuple(float(v) for v in incidence_range)
        except (TypeError, ValueError):
            raise ValueError("incidence_range must contain two numeric values")
        if len(incidence_range_tuple) != 2:
            raise ValueError("incidence_range must be a 2-element tuple (near, far)")

        reference_incidence_angle = 0.5 * (
            incidence_range_tuple[0] + incidence_range_tuple[1]
        )

        dem_quality_flags = self.quality_metrics.get('dem_quality_flags')
        if dem_quality_flags is not None and not isinstance(dem_quality_flags, dict):
            raise ValueError("dem_quality_flags quality metric must be a dictionary if provided")

        dem_void_fill = self.quality_metrics.get('dem_void_fill_method')
        if dem_void_fill is not None and not isinstance(dem_void_fill, str):
            raise ValueError("dem_void_fill_method quality metric must be a string if provided")

        dem_metadata = DEMMetadata(
            source=dem_source,
            version=dem_version,
            resolution_meters=float(dem_resolution),
            vertical_datum="EGM96",
            horizontal_datum="WGS84",
            coverage_area=dem_coverage_tuple,
            void_fill_method=dem_void_fill,
            smoothing_applied=bool(self.quality_metrics.get('dem_smoothing_applied', False)),
            quality_flags=dem_quality_flags,
        )

        rtc_processing = RTCProcessingMetadata(
            rtc_method=rtc_method,
            correction_formula="γ⁰ = σ⁰ / cos(θ_local)",
            local_incidence_method="DEM_normals",
            cosine_clip_threshold=float(cosine_threshold),
            masking_applied=bool(masking_enabled),
            shadow_layover_detection=bool(
                shadow_pixels_pct is not None or layover_pixels_pct is not None
            ),
            reference_incidence_angle=reference_incidence_angle,
            processing_mode="linear",
        )

        calibration_metadata = CalibrationMetadata(
            calibration_type=calibration_type,
            calibration_equation="σ⁰ = β⁰ · sin(θ_incidence)",
            lut_source=lut_source,
            noise_removal_applied=bool(noise_applied),
            noise_source="annotation_derived",
            antenna_pattern_correction=True,
            range_spreading_correction=True,
        )
        
        # Build processing chain metadata with provenance (PROV-1/PROV-2)
        provenance_dict = _get_provenance()
        processing_chain = ProcessingChainMetadata(
            processing_steps=self.processing_steps,
            software_version=self.software_version,
            processing_timestamp=datetime.now(timezone.utc).isoformat(),
            input_product=self.product_id,
            orbit_files_applied=self.orbit_files,
            total_processing_time=total_time,
            provenance=provenance_dict,
        )
        
        # Build quality metadata
        quality_flags: List[str] = []

        quality_metadata = QualityMetadata(
            valid_pixel_percentage=float(self.quality_metrics['valid_pixel_percentage']),
            shadow_pixel_percentage=shadow_pixels_pct,
            layover_pixel_percentage=layover_pixels_pct,
            steep_terrain_percentage=steep_terrain_pct,
            backscatter_statistics=backscatter_stats,
            incidence_angle_range=incidence_range_tuple,
            quality_flags=quality_flags
        )
        
        # Build geospatial metadata
        geospatial_metadata = GeospatialMetadata(
            coordinate_system=coordinate_system,
            pixel_spacing=pixel_spacing_tuple,
            scene_center_lat=scene_center_tuple[1],
            scene_center_lon=scene_center_tuple[0],
            scene_extent=scene_extent_tuple,
            acquisition_geometry={
                'look_direction': 'right',
                'polarization': self.polarization,
                'incidence_angle_range': incidence_range_tuple
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
    
    # Record processing artefacts required for metadata
    builder.set_software_version("SARdine v0.2.0")
    builder.add_orbit_file("S1A_OPER_AUX_POEORB_OPOD_20240720T120000_V20240719T235942_20240721T015942.EOF")
    builder.set_quality_metrics(
        {
            "valid_pixel_percentage": 95.2,
            "backscatter_statistics": {"min": -30.5, "max": 15.2, "mean": -12.8, "std": 4.3},
            "shadow_pixel_percentage": 2.4,
            "layover_pixel_percentage": 1.2,
            "steep_terrain_percentage": 0.6,
        }
    )

    # Build comprehensive metadata
    metadata = builder.build_metadata(
        dem_source="SRTM_30m",
        dem_resolution=30.0,
        dem_version="v3.0",
        dem_coverage=(-122.5, 37.5, -122.0, 38.0),
        rtc_method="gamma_flattening_with_masking",
        cosine_threshold=0.1,
        masking_enabled=True,
        calibration_type="sigma0",
        lut_source="calibration.xml",
        noise_applied=True,
        incidence_range=(29.1, 45.9),
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
