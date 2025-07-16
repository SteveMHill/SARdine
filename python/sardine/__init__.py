"""
SARdine: A Fast, Modular Sentinel-1 Backscatter Processor

A modern, open-source alternative to ESA SNAP and GAMMA for processing 
Sentinel-1 SLC data into calibrated, terrain-corrected backscatter products.
"""

from . import _core
from ._core import (
    SlcReader, 
    Polarization, 
    Metadata, 
    test_srtm_download, 
    apply_speckle_filter, 
    estimate_num_looks
)

__version__ = "0.1.0"
__author__ = "Steven Hill and contributors"
__license__ = "MIT"

__all__ = [
    "SlcReader",
    "Polarization", 
    "Metadata",
    "test_srtm_download",
    "apply_speckle_filter",
    "estimate_num_looks",
    "process_slc",
    "get_product_info",
]

def process_slc(input_path, output_path=None, **kwargs):
    """
    High-level function to process Sentinel-1 SLC data.
    
    Parameters
    ----------
    input_path : str
        Path to Sentinel-1 SLC ZIP file
    output_path : str, optional
        Output directory for processed products
    **kwargs
        Additional processing parameters
        
    Returns
    -------
    dict
        Processing results and metadata
    """
    # This will be implemented as we add more processing modules
    raise NotImplementedError("Full processing pipeline not yet implemented")

def get_product_info(input_path):
    """
    Extract and display information about a Sentinel-1 product.
    
    Parameters
    ----------
    input_path : str
        Path to Sentinel-1 SLC ZIP file
        
    Returns
    -------
    dict
        Product information and metadata
    """
    reader = SlcReader(input_path)
    files = reader.list_files()
    
    # Find available polarizations
    polarizations = []
    if any("vv" in f.lower() for f in files):
        polarizations.append("VV")
    if any("vh" in f.lower() for f in files):
        polarizations.append("VH")
    if any("hv" in f.lower() for f in files):
        polarizations.append("HV")
    if any("hh" in f.lower() for f in files):
        polarizations.append("HH")
    
    result = {
        "input_path": input_path,
        "total_files": len(files),
        "polarizations": polarizations,
        "metadata": {}
    }
    
    # Get metadata for first available polarization
    if polarizations:
        try:
            metadata = reader.get_metadata(polarizations[0])
            result["metadata"] = {
                "product_id": metadata.product_id,
                "mission": metadata.mission,
                "platform": metadata.platform,
                "start_time": metadata.start_time,
                "stop_time": metadata.stop_time,
                "acquisition_mode": metadata.acquisition_mode,
                "pixel_spacing": metadata.pixel_spacing,
                "bounding_box": metadata.bounding_box,
            }
        except Exception as e:
            result["metadata_error"] = str(e)
    
    return result
