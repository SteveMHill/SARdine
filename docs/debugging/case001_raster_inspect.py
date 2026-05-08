#!/usr/bin/env python3
"""
Temporary artifact forensics script for Case 001.
Reads the processed VV backscatter raster and reports diagnostics.
Does NOT modify any data or pipeline code.

Usage: python3 docs/debugging/case001_raster_inspect.py
"""

import sys
import numpy as np
import rasterio

RASTER_PATH = "/home/datacube/dev/SARdine/docs/debugging/backscatter_VV_final.tif"


def main():
    print("=" * 72)
    print("CASE 001 ARTIFACT FORENSICS — RASTER INSPECTION")
    print("=" * 72)

    with rasterio.open(RASTER_PATH) as ds:
        # ── Section 1: Basic metadata ──
        print("\n── 1. RASTER METADATA ──")
        print(f"Path:       {RASTER_PATH}")
        print(f"Driver:     {ds.driver}")
        print(f"Shape:      {ds.height} rows × {ds.width} cols  (bands={ds.count})")
        print(f"Dtype:      {ds.dtypes[0]}")
        print(f"CRS:        {ds.crs}")
        print(f"Transform:  {ds.transform}")
        print(f"Bounds:     {ds.bounds}")
        print(f"Nodata:     {ds.nodata}")
        total_pixels = ds.height * ds.width
        print(f"Total px:   {total_pixels:,}")

        # ── Section 2: Read data and compute valid/invalid ──
        print("\n── 2. VALID / INVALID PIXEL COUNTS ──")
        data = ds.read(1)  # Band 1
        nodata_val = ds.nodata

        # Invalid = NaN OR matches nodata sentinel
        is_nan = np.isnan(data)
        if nodata_val is not None and not np.isnan(nodata_val):
            is_nodata = (data == nodata_val)
            invalid_mask = is_nan | is_nodata
        else:
            # nodata is NaN or not set — NaN is the sole invalid marker
            invalid_mask = is_nan

        valid_mask = ~invalid_mask
        valid_count = int(np.sum(valid_mask))
        invalid_count = int(np.sum(invalid_mask))
        print(f"Valid px:    {valid_count:,}  ({100*valid_count/total_pixels:.2f}%)")
        print(f"Invalid px: {invalid_count:,}  ({100*invalid_count/total_pixels:.2f}%)")

        # ── Section 3: Value statistics on valid pixels ──
        print("\n── 3. VALUE STATISTICS (valid pixels only) ──")
        valid_vals = data[valid_mask]
        if len(valid_vals) > 0:
            print(f"Min:        {float(np.min(valid_vals)):.6f}")
            print(f"Max:        {float(np.max(valid_vals)):.6f}")
            print(f"Mean:       {float(np.mean(valid_vals)):.6f}")
            print(f"Median:     {float(np.median(valid_vals)):.6f}")
            print(f"Std:        {float(np.std(valid_vals)):.6f}")
            # Check for negative values (expected in dB domain)
            neg_count = int(np.sum(valid_vals < 0))
            pos_count = int(np.sum(valid_vals > 0))
            zero_count = int(np.sum(valid_vals == 0))
            print(f"Negative:   {neg_count:,}  Positive: {pos_count:,}  Zero: {zero_count:,}")
        else:
            print("(no valid pixels)")

        # ── Section 4: Row-level invalid structure ──
        print("\n── 4. ROW-LEVEL INVALID STRUCTURE ──")
        invalid_per_row = np.sum(invalid_mask, axis=1)  # shape: (height,)
        fully_invalid_rows = int(np.sum(invalid_per_row == ds.width))
        fully_valid_rows = int(np.sum(invalid_per_row == 0))
        partial_rows = ds.height - fully_invalid_rows - fully_valid_rows
        print(f"Fully invalid rows:  {fully_invalid_rows}")
        print(f"Fully valid rows:    {fully_valid_rows}")
        print(f"Partially valid:     {partial_rows}")

        # Find contiguous bands of fully-invalid rows
        is_full_invalid_row = (invalid_per_row == ds.width)
        bands = _find_contiguous_runs(is_full_invalid_row)
        if bands:
            print(f"\nContiguous fully-invalid row bands ({len(bands)} found):")
            for start, length in bands[:20]:  # Show first 20
                print(f"  rows {start}–{start+length-1}  (length={length})")
            if len(bands) > 20:
                print(f"  ... and {len(bands)-20} more bands")
        else:
            print("\nNo contiguous fully-invalid row bands found.")

        # ── Section 5: Column-level invalid structure ──
        print("\n── 5. COLUMN-LEVEL INVALID STRUCTURE ──")
        invalid_per_col = np.sum(invalid_mask, axis=0)  # shape: (width,)
        fully_invalid_cols = int(np.sum(invalid_per_col == ds.height))
        fully_valid_cols = int(np.sum(invalid_per_col == 0))
        partial_cols = ds.width - fully_invalid_cols - fully_valid_cols
        print(f"Fully invalid cols:  {fully_invalid_cols}")
        print(f"Fully valid cols:    {fully_valid_cols}")
        print(f"Partially valid:     {partial_cols}")

        col_bands = _find_contiguous_runs(invalid_per_col == ds.height)
        if col_bands:
            print(f"\nContiguous fully-invalid column bands ({len(col_bands)} found):")
            for start, length in col_bands[:20]:
                print(f"  cols {start}–{start+length-1}  (length={length})")
        else:
            print("\nNo contiguous fully-invalid column bands found.")

        # ── Section 6: Invalid fraction per row — profile ──
        print("\n── 6. INVALID FRACTION PER ROW — PROFILE (sampled) ──")
        invalid_frac_per_row = invalid_per_row / ds.width
        # Sample every Nth row to show the profile compactly
        step = max(1, ds.height // 40)
        print(f"  {'Row':>6s}  {'Invalid%':>8s}  Bar")
        for r in range(0, ds.height, step):
            frac = invalid_frac_per_row[r]
            bar = "#" * int(frac * 50)
            print(f"  {r:6d}  {frac*100:7.1f}%  {bar}")

        # ── Section 7: Row-mean profile for duplication detection ──
        print("\n── 7. ROW-MEAN PROFILE — DUPLICATION DETECTION ──")
        # Compute mean of valid pixels per row (NaN where fully invalid)
        with np.errstate(invalid='ignore'):
            row_means = np.where(
                invalid_per_row < ds.width,
                np.nanmean(np.where(valid_mask, data, np.nan), axis=1),
                np.nan
            )

        # Look for pairs of rows with nearly identical mean
        # Focus on close-by rows (potential burst-period duplicates)
        # Typical burst height in IW mode post-multilook: ~300-1500 rows
        print("Checking for repeating row-mean patterns at burst-like periods...")
        valid_row_means = row_means[np.isfinite(row_means)]
        if len(valid_row_means) < 100:
            print("  Too few valid rows for duplication analysis.")
        else:
            # Autocorrelation at candidate lags
            centered = valid_row_means - np.mean(valid_row_means)
            norm = np.sum(centered**2)
            if norm > 0:
                candidate_lags = [100, 200, 300, 400, 500, 600, 700, 800,
                                  900, 1000, 1200, 1500, 2000, 2500, 3000]
                print(f"  {'Lag':>6s}  {'AutoCorr':>9s}")
                for lag in candidate_lags:
                    if lag < len(centered):
                        ac = np.sum(centered[:-lag] * centered[lag:]) / norm
                        print(f"  {lag:6d}  {ac:9.4f}")

        # ── Section 8: Direct row-pair comparison ──
        print("\n── 8. DIRECT ROW-PAIR DUPLICATION CHECK ──")
        # Pick a central valid row, compare it to rows at various offsets
        valid_row_indices = np.where(invalid_per_row < ds.width * 0.5)[0]
        if len(valid_row_indices) > 0:
            mid_idx = len(valid_row_indices) // 2
            ref_row_idx = valid_row_indices[mid_idx]
            ref_row = data[ref_row_idx, :].astype(np.float64)
            ref_valid = valid_mask[ref_row_idx, :]

            print(f"Reference row: {ref_row_idx}")
            offsets_to_check = [1, 2, 5, 10, 50, 100, 200, 300, 500, 800, 1000, 1500, 2000]
            print(f"  {'Offset':>7s}  {'CompRow':>7s}  {'SharedPx':>9s}  {'MAE':>10s}  {'MaxDiff':>10s}  {'Corr':>7s}")
            for offset in offsets_to_check:
                comp_row_idx = ref_row_idx + offset
                if comp_row_idx >= ds.height:
                    continue
                comp_row = data[comp_row_idx, :].astype(np.float64)
                comp_valid = valid_mask[comp_row_idx, :]
                shared = ref_valid & comp_valid
                n_shared = int(np.sum(shared))
                if n_shared < 100:
                    continue
                diff = np.abs(ref_row[shared] - comp_row[shared])
                mae = float(np.mean(diff))
                max_diff = float(np.max(diff))
                # Pearson correlation
                a = ref_row[shared]
                b = comp_row[shared]
                a_c = a - np.mean(a)
                b_c = b - np.mean(b)
                denom = np.sqrt(np.sum(a_c**2) * np.sum(b_c**2))
                corr = float(np.sum(a_c * b_c) / denom) if denom > 0 else 0.0
                print(f"  {offset:7d}  {comp_row_idx:7d}  {n_shared:9d}  {mae:10.6f}  {max_diff:10.6f}  {corr:7.4f}")

        # ── Section 9: Edge vs interior invalid analysis ──
        print("\n── 9. EDGE VS INTERIOR INVALID ANALYSIS ──")
        # Check if invalids are concentrated at edges (expected for rotated
        # SAR footprint in projected CRS) or interior (structural defect)
        edge_width = min(200, ds.width // 10, ds.height // 10)
        top_edge = invalid_mask[:edge_width, :]
        bottom_edge = invalid_mask[-edge_width:, :]
        left_edge = invalid_mask[:, :edge_width]
        right_edge = invalid_mask[:, -edge_width:]
        interior = invalid_mask[edge_width:-edge_width, edge_width:-edge_width]

        top_frac = np.mean(top_edge)
        bot_frac = np.mean(bottom_edge)
        left_frac = np.mean(left_edge)
        right_frac = np.mean(right_edge)
        int_frac = np.mean(interior)
        print(f"Top {edge_width} rows:     {top_frac*100:.1f}% invalid")
        print(f"Bottom {edge_width} rows:  {bot_frac*100:.1f}% invalid")
        print(f"Left {edge_width} cols:    {left_frac*100:.1f}% invalid")
        print(f"Right {edge_width} cols:   {right_frac*100:.1f}% invalid")
        print(f"Interior:           {int_frac*100:.1f}% invalid")

        if int_frac > 0.01:
            print("\n⚠  Interior invalid fraction > 1% — possible structural defect")
            # Find interior invalid clusters
            interior_inv = invalid_mask[edge_width:-edge_width, edge_width:-edge_width]
            inv_rows_int = np.sum(interior_inv, axis=1)
            high_inv_rows = np.where(inv_rows_int > interior_inv.shape[1] * 0.5)[0]
            if len(high_inv_rows) > 0:
                bands_int = _find_contiguous_runs(inv_rows_int > interior_inv.shape[1] * 0.5)
                print(f"  Interior row bands with >50% invalid ({len(bands_int)} found):")
                for start, length in bands_int[:15]:
                    actual_row = start + edge_width
                    print(f"    rows {actual_row}–{actual_row+length-1}  (length={length})")

        # ── Section 10: Summary ──
        print("\n" + "=" * 72)
        print("CASE 001 ARTIFACT FORENSICS — SUMMARY")
        print("=" * 72)
        print(f"Shape:              {ds.height} × {ds.width}")
        print(f"Valid coverage:     {100*valid_count/total_pixels:.2f}%")
        print(f"Invalid structure:  ", end="")
        if fully_invalid_rows > 0 and len(bands) > 0:
            print(f"{len(bands)} contiguous invalid row-bands detected")
        elif fully_invalid_cols > 0:
            print(f"{fully_invalid_cols} fully invalid columns detected")
        else:
            print("no full row/column bands; scattered or edge-only")
        print(f"Edge vs interior:   edge={max(top_frac,bot_frac,left_frac,right_frac)*100:.0f}%, interior={int_frac*100:.1f}%")

    print("\nDone.")


def _find_contiguous_runs(bool_array):
    """Find contiguous True runs. Returns list of (start_index, length)."""
    runs = []
    in_run = False
    start = 0
    for i, v in enumerate(bool_array):
        if v and not in_run:
            in_run = True
            start = i
        elif not v and in_run:
            runs.append((start, i - start))
            in_run = False
    if in_run:
        runs.append((start, len(bool_array) - start))
    return runs


if __name__ == "__main__":
    main()
