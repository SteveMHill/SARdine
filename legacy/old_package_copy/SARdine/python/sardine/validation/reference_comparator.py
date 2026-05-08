"""
Reference Product Comparator for SARdine

Compares SARdine outputs against reference products from:
- ASF (Alaska Satellite Facility) RTC products
- SNAP (ESA Sentinel Application Platform)
- Other calibrated reference sources

This enables validation that SARdine produces scientifically correct results.
"""

import logging
import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional, Tuple, Union
from pathlib import Path

try:
    import rasterio
    from rasterio.warp import reproject, Resampling
    from rasterio.crs import CRS
    HAS_RASTERIO = True
except ImportError:
    HAS_RASTERIO = False

logger = logging.getLogger(__name__)


@dataclass
class ComparisonResult:
    """Result of comparing two products."""
    metric: str
    value: float
    passed: bool
    threshold: float
    details: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ReferenceComparisonReport:
    """Complete comparison report."""
    sardine_path: str
    reference_path: str
    reference_source: str  # "ASF", "SNAP", "other"
    timestamp: str
    results: List[ComparisonResult] = field(default_factory=list)
    spatial_overlap_percent: float = 0.0
    common_valid_pixels: int = 0
    
    @property
    def passed(self) -> bool:
        """Check if all comparisons passed."""
        return all(r.passed for r in self.results)
    
    @property
    def mean_absolute_error_db(self) -> Optional[float]:
        """Get MAE in dB if available."""
        for r in self.results:
            if r.metric == "mean_absolute_error_db":
                return r.value
        return None


class ReferenceComparator:
    """
    Compare SARdine outputs against reference products.
    
    Supports comparison with:
    - ASF RTC products (gamma0 terrain-corrected)
    - SNAP processed outputs
    - Any GeoTIFF with compatible coverage
    
    Metrics computed:
    - Mean Absolute Error (MAE) in dB
    - Root Mean Square Error (RMSE) in dB
    - Correlation coefficient
    - Bias (systematic offset)
    - Spatial statistics
    """
    
    # Thresholds based on ESA CCI Land Cover validation spec
    DEFAULT_THRESHOLDS = {
        "max_mae_db": 1.0,           # Maximum allowed mean absolute error
        "max_rmse_db": 1.5,          # Maximum allowed RMSE
        "min_correlation": 0.95,      # Minimum required correlation
        "max_bias_db": 0.5,          # Maximum allowed systematic bias
        "max_std_diff_db": 2.0,      # Maximum std difference
        "min_overlap_percent": 50.0,  # Minimum spatial overlap
    }
    
    def __init__(self, thresholds: Optional[Dict[str, float]] = None):
        """
        Initialize comparator.
        
        Args:
            thresholds: Custom thresholds (overrides defaults)
        """
        self.thresholds = self.DEFAULT_THRESHOLDS.copy()
        if thresholds:
            self.thresholds.update(thresholds)
    
    def compare(
        self,
        sardine_path: Union[str, Path],
        reference_path: Union[str, Path],
        reference_source: str = "unknown",
        band_sardine: int = 1,
        band_reference: int = 1,
    ) -> ReferenceComparisonReport:
        """
        Compare SARdine output with reference product.
        
        Args:
            sardine_path: Path to SARdine output GeoTIFF
            reference_path: Path to reference product GeoTIFF
            reference_source: Source of reference ("ASF", "SNAP", etc.)
            band_sardine: Band number in SARdine product (1-indexed)
            band_reference: Band number in reference product (1-indexed)
            
        Returns:
            ReferenceComparisonReport with all comparison metrics
        """
        from datetime import datetime
        
        if not HAS_RASTERIO:
            raise ImportError("rasterio is required for reference comparison")
        
        report = ReferenceComparisonReport(
            sardine_path=str(sardine_path),
            reference_path=str(reference_path),
            reference_source=reference_source,
            timestamp=datetime.now().isoformat(),
        )
        
        try:
            # Load products and align them
            sardine_data, reference_data, overlap_percent = self._load_and_align(
                sardine_path, reference_path, band_sardine, band_reference
            )
            
            report.spatial_overlap_percent = overlap_percent
            
            # Check minimum overlap
            if overlap_percent < self.thresholds["min_overlap_percent"]:
                report.results.append(ComparisonResult(
                    metric="spatial_overlap",
                    value=overlap_percent,
                    passed=False,
                    threshold=self.thresholds["min_overlap_percent"],
                    details={"message": "Insufficient spatial overlap"}
                ))
                return report
            
            # Convert to dB if needed
            sardine_db = self._to_db(sardine_data)
            reference_db = self._to_db(reference_data)
            
            # Find common valid pixels
            valid_mask = (
                np.isfinite(sardine_db) & 
                np.isfinite(reference_db) &
                (sardine_db > -50) & (sardine_db < 30) &
                (reference_db > -50) & (reference_db < 30)
            )
            
            report.common_valid_pixels = int(np.sum(valid_mask))
            
            if report.common_valid_pixels < 1000:
                report.results.append(ComparisonResult(
                    metric="valid_pixels",
                    value=float(report.common_valid_pixels),
                    passed=False,
                    threshold=1000.0,
                    details={"message": "Too few valid pixels for comparison"}
                ))
                return report
            
            # Extract valid pixels
            s_valid = sardine_db[valid_mask]
            r_valid = reference_db[valid_mask]
            
            # Compute comparison metrics
            report.results.extend(self._compute_metrics(s_valid, r_valid))
            
            # Add spatial overlap result
            report.results.append(ComparisonResult(
                metric="spatial_overlap",
                value=overlap_percent,
                passed=overlap_percent >= self.thresholds["min_overlap_percent"],
                threshold=self.thresholds["min_overlap_percent"],
            ))
            
        except Exception as e:
            logger.error(f"Comparison failed: {e}")
            report.results.append(ComparisonResult(
                metric="comparison_error",
                value=0.0,
                passed=False,
                threshold=0.0,
                details={"error": str(e)}
            ))
        
        return report
    
    def _load_and_align(
        self,
        sardine_path: Union[str, Path],
        reference_path: Union[str, Path],
        band_sardine: int,
        band_reference: int,
    ) -> Tuple[np.ndarray, np.ndarray, float]:
        """
        Load and spatially align two products.
        
        Returns:
            Tuple of (sardine_data, reference_data, overlap_percent)
        """
        with rasterio.open(sardine_path) as src_s, rasterio.open(reference_path) as src_r:
            # Get transforms and shapes
            transform_s = src_s.transform
            transform_r = src_r.transform
            
            crs_s = src_s.crs
            crs_r = src_r.crs
            
            # Calculate overlapping bounds
            bounds_s = src_s.bounds
            bounds_r = src_r.bounds
            
            # If different CRS, transform reference bounds to sardine CRS
            if crs_s != crs_r:
                from rasterio.warp import transform_bounds
                bounds_r = transform_bounds(crs_r, crs_s, *bounds_r)
            
            # Find intersection
            overlap_left = max(bounds_s.left, bounds_r.left)
            overlap_bottom = max(bounds_s.bottom, bounds_r.bottom)
            overlap_right = min(bounds_s.right, bounds_r.right)
            overlap_top = min(bounds_s.top, bounds_r.top)
            
            if overlap_left >= overlap_right or overlap_bottom >= overlap_top:
                return np.array([]), np.array([]), 0.0
            
            # Calculate overlap percentage
            sardine_area = (bounds_s.right - bounds_s.left) * (bounds_s.top - bounds_s.bottom)
            overlap_area = (overlap_right - overlap_left) * (overlap_top - overlap_bottom)
            overlap_percent = overlap_area / sardine_area * 100.0
            
            # Read sardine data for overlap region
            window_s = rasterio.windows.from_bounds(
                overlap_left, overlap_bottom, overlap_right, overlap_top,
                transform_s
            )
            sardine_data = src_s.read(band_sardine, window=window_s)
            
            # Get the actual transform for the window
            window_transform = rasterio.windows.transform(window_s, transform_s)
            
            # Read and reproject reference data to match
            reference_data = np.zeros_like(sardine_data)
            
            reproject(
                source=rasterio.band(src_r, band_reference),
                destination=reference_data,
                src_transform=transform_r,
                src_crs=crs_r,
                dst_transform=window_transform,
                dst_crs=crs_s,
                resampling=Resampling.bilinear
            )
            
            return sardine_data, reference_data, overlap_percent
    
    def _to_db(self, data: np.ndarray) -> np.ndarray:
        """
        Convert data to dB if it appears to be linear.
        
        Heuristic: if data is mostly positive and < 100, it's probably linear.
        """
        finite = data[np.isfinite(data)]
        if len(finite) == 0:
            return data
        
        data_min = np.min(finite)
        data_max = np.max(finite)
        
        # Appears to be linear power data
        if data_min >= 0 and data_max < 100:
            # Convert to dB, avoiding log(0)
            with np.errstate(divide='ignore', invalid='ignore'):
                db_data = 10.0 * np.log10(np.maximum(data, 1e-10))
            return db_data
        
        # Already in dB (has negative values or typical dB range)
        return data
    
    def _compute_metrics(
        self, 
        sardine: np.ndarray, 
        reference: np.ndarray
    ) -> List[ComparisonResult]:
        """Compute all comparison metrics."""
        results = []
        
        # 1. Mean Absolute Error
        diff = sardine - reference
        mae = float(np.mean(np.abs(diff)))
        results.append(ComparisonResult(
            metric="mean_absolute_error_db",
            value=mae,
            passed=mae <= self.thresholds["max_mae_db"],
            threshold=self.thresholds["max_mae_db"],
            details={
                "interpretation": "Average absolute difference in dB",
            }
        ))
        
        # 2. Root Mean Square Error
        rmse = float(np.sqrt(np.mean(diff ** 2)))
        results.append(ComparisonResult(
            metric="rmse_db",
            value=rmse,
            passed=rmse <= self.thresholds["max_rmse_db"],
            threshold=self.thresholds["max_rmse_db"],
            details={
                "interpretation": "RMS difference emphasizing large errors",
            }
        ))
        
        # 3. Bias (mean difference)
        bias = float(np.mean(diff))
        results.append(ComparisonResult(
            metric="bias_db",
            value=bias,
            passed=abs(bias) <= self.thresholds["max_bias_db"],
            threshold=self.thresholds["max_bias_db"],
            details={
                "interpretation": "Systematic offset (positive = SARdine brighter)",
            }
        ))
        
        # 4. Correlation coefficient
        correlation = float(np.corrcoef(sardine.flatten(), reference.flatten())[0, 1])
        results.append(ComparisonResult(
            metric="correlation",
            value=correlation,
            passed=correlation >= self.thresholds["min_correlation"],
            threshold=self.thresholds["min_correlation"],
            details={
                "interpretation": "Spatial pattern agreement",
            }
        ))
        
        # 5. Standard deviation difference
        std_sardine = float(np.std(sardine))
        std_reference = float(np.std(reference))
        std_diff = abs(std_sardine - std_reference)
        results.append(ComparisonResult(
            metric="std_difference_db",
            value=std_diff,
            passed=std_diff <= self.thresholds["max_std_diff_db"],
            threshold=self.thresholds["max_std_diff_db"],
            details={
                "std_sardine": std_sardine,
                "std_reference": std_reference,
                "interpretation": "Difference in dynamic range/contrast",
            }
        ))
        
        # 6. Percentile comparison
        percentiles = [5, 25, 50, 75, 95]
        sardine_percentiles = np.percentile(sardine, percentiles)
        reference_percentiles = np.percentile(reference, percentiles)
        percentile_diffs = sardine_percentiles - reference_percentiles
        
        results.append(ComparisonResult(
            metric="percentile_analysis",
            value=float(np.max(np.abs(percentile_diffs))),
            passed=np.all(np.abs(percentile_diffs) < 2.0),  # 2 dB tolerance
            threshold=2.0,
            details={
                "percentiles": percentiles,
                "sardine_values": sardine_percentiles.tolist(),
                "reference_values": reference_percentiles.tolist(),
                "differences": percentile_diffs.tolist(),
            }
        ))
        
        return results
    
    def compare_asf_rtc(
        self,
        sardine_path: Union[str, Path],
        asf_rtc_path: Union[str, Path],
    ) -> ReferenceComparisonReport:
        """
        Compare with ASF RTC (Radiometric Terrain Corrected) product.
        
        ASF RTC products are gamma0 with 30m resolution in UTM.
        
        Args:
            sardine_path: Path to SARdine output
            asf_rtc_path: Path to ASF RTC product (VV or VH band)
            
        Returns:
            Comparison report
        """
        return self.compare(
            sardine_path,
            asf_rtc_path,
            reference_source="ASF_RTC",
        )
    
    def compare_snap(
        self,
        sardine_path: Union[str, Path],
        snap_path: Union[str, Path],
    ) -> ReferenceComparisonReport:
        """
        Compare with SNAP processed product.
        
        Args:
            sardine_path: Path to SARdine output
            snap_path: Path to SNAP output GeoTIFF
            
        Returns:
            Comparison report
        """
        return self.compare(
            sardine_path,
            snap_path,
            reference_source="SNAP",
        )
    
    def print_report(self, report: ReferenceComparisonReport):
        """Print formatted comparison report."""
        print("\n" + "=" * 80)
        print("REFERENCE COMPARISON REPORT")
        print("=" * 80)
        print(f"SARdine product: {report.sardine_path}")
        print(f"Reference product: {report.reference_path}")
        print(f"Reference source: {report.reference_source}")
        print(f"Timestamp: {report.timestamp}")
        print()
        print(f"Spatial overlap: {report.spatial_overlap_percent:.1f}%")
        print(f"Common valid pixels: {report.common_valid_pixels:,}")
        print()
        
        print("Overall: " + ("✅ PASSED" if report.passed else "❌ FAILED"))
        print()
        
        print("METRICS:")
        print("-" * 60)
        for r in report.results:
            status = "✅" if r.passed else "❌"
            
            if r.metric == "bias_db":
                # Show bias with sign
                print(f"{status} {r.metric}: {r.value:+.3f} dB (threshold: ±{r.threshold} dB)")
            elif "db" in r.metric.lower():
                print(f"{status} {r.metric}: {r.value:.3f} dB (threshold: {r.threshold} dB)")
            elif r.metric == "correlation":
                print(f"{status} {r.metric}: {r.value:.4f} (threshold: ≥{r.threshold})")
            elif r.metric == "spatial_overlap":
                print(f"{status} {r.metric}: {r.value:.1f}% (threshold: ≥{r.threshold}%)")
            else:
                print(f"{status} {r.metric}: {r.value:.3f} (threshold: {r.threshold})")
            
            if "interpretation" in r.details:
                print(f"     ℹ️  {r.details['interpretation']}")
        
        print("=" * 80)
    
    def generate_difference_map(
        self,
        sardine_path: Union[str, Path],
        reference_path: Union[str, Path],
        output_path: Union[str, Path],
    ) -> bool:
        """
        Generate a difference map GeoTIFF for visual inspection.
        
        Args:
            sardine_path: Path to SARdine output
            reference_path: Path to reference product
            output_path: Path for output difference map
            
        Returns:
            True if successful
        """
        try:
            sardine_data, reference_data, overlap = self._load_and_align(
                sardine_path, reference_path, 1, 1
            )
            
            if overlap < 10:
                logger.error("Insufficient overlap for difference map")
                return False
            
            # Convert to dB
            sardine_db = self._to_db(sardine_data)
            reference_db = self._to_db(reference_data)
            
            # Calculate difference
            diff = sardine_db - reference_db
            
            # Write output
            with rasterio.open(sardine_path) as src:
                profile = src.profile.copy()
                profile.update(dtype='float32', count=1)
                
                with rasterio.open(output_path, 'w', **profile) as dst:
                    dst.write(diff.astype(np.float32), 1)
            
            logger.info(f"Difference map written to {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"Failed to generate difference map: {e}")
            return False
