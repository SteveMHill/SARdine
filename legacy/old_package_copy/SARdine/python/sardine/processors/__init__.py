"""
SARdine Processors Package

Contains specialized processors for different SAR processing workflows.
"""

from .backscatter import BackscatterProcessor
from .base import BaseProcessor

__all__ = ['BackscatterProcessor', 'BaseProcessor']
