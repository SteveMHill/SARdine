"""
SARdine: High-Performance SAR Processing Library

A modern SAR data processing library for Sentinel-1 data, implemented in Rust 
with Python bindings. Provides complete processing pipeline from SLC to 
analysis-ready backscatter products.
"""

__version__ = "0.2.0"
__author__ = "SARdine Contributors"

# Import all functions from the core module
from ._core import *

# Feature availability flags
_TERRAIN_CORRECTION_AVAILABLE = True
_TERRAIN_FLATTENING_AVAILABLE = True
