#!/usr/bin/env python3
"""
Case 001 deep artifact analysis — Part 3.
Focused on interior invalid pixel structure:
  - In rows 2500–7500 (well-covered mid-scene), ~14% are still invalid
  - Are these isolated pixels, small holes, or structured bands?
  - Where along each row do the interior invalids appear?
"""

import numpy as np
import rasterio

RASTER_PATH = "/home/datacube/dev/SARdine/docs/debugging/backscatter_VV_final.tif"


def main():
    with rasterio.open(RASTER_PATH) as ds:
        data = ds.read(1)
        invalid = np.isnan(data)
        valid = ~invalid
        H, W = data.shape

        # ═══════ A: Interior invalid distribution ═══════
        print("=" * 72)
        print("A. INTERIOR INVALID PIXEL COLUMN DISTRIBUTION (rows 3000–7500)")
        print("=" * 72)
        # For rows in the well-covered middle, where along the row are invalids?
        interior_rows = range(3000, 7500)
        inv_col_histogram = np.zeros(W, dtype=np.int64)
        for r in interior_rows:
            inv_col_histogram += invalid[r, :].astype(np.int64)

        n_rows_sampled = len(interior_rows)
        print(f"Rows sampled: {n_rows_sampled}")
        print(f"\nColumn-binned invalid fraction (100-col bins):")
        bin_size = 100
        print(f"  {'ColRange':>12s}  {'InvFrac':>8s}  {'Bar'}")
        for c_start in range(0, W, bin_size):
            c_end = min(c_start + bin_size, W)
            frac = np.sum(inv_col_histogram[c_start:c_end]) / (n_rows_sampled * (c_end - c_start))
            bar = "#" * int(frac * 50)
            if frac > 0.01:
                print(f"  {c_start:5d}–{c_end:5d}  {frac*100:7.1f}%  {bar}")

        # ═══════ B: Invalid run-length analysis in interior rows ═══════
        print("\n" + "=" * 72)
        print("B. INVALID RUN-LENGTH ANALYSIS IN INTERIOR ROWS")
        print("=" * 72)
        # For each row in range, find contiguous invalid runs
        all_runs = []  # (row, col_start, length)
        for r in range(3000, 7500, 10):  # Sample every 10th row
            inv_row = invalid[r, :]
            # Find first/last valid to exclude edge padding
            v_cols = np.where(~inv_row)[0]
            if len(v_cols) < 100:
                continue
            first_valid = v_cols[0]
            last_valid = v_cols[-1]
            # Find invalid runs within the valid range
            in_run = False
            run_start = 0
            for c in range(first_valid, last_valid + 1):
                if inv_row[c] and not in_run:
                    in_run = True
                    run_start = c
                elif not inv_row[c] and in_run:
                    all_runs.append((r, run_start, c - run_start))
                    in_run = False
            if in_run:
                all_runs.append((r, run_start, last_valid + 1 - run_start))

        if all_runs:
            lengths = [l for _, _, l in all_runs]
            print(f"Total interior invalid runs found: {len(all_runs)}")
            print(f"Min run length: {min(lengths)}")
            print(f"Max run length: {max(lengths)}")
            print(f"Mean run length: {np.mean(lengths):.1f}")
            print(f"Median run length: {np.median(lengths):.1f}")

            # Histogram of run lengths
            print(f"\nRun length distribution:")
            bins = [1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 5000, 20000]
            for i in range(len(bins) - 1):
                count = sum(1 for l in lengths if bins[i] <= l < bins[i + 1])
                if count > 0:
                    print(f"  {bins[i]:5d}–{bins[i+1]-1:5d}: {count:6d} runs")
            count_large = sum(1 for l in lengths if l >= bins[-1])
            if count_large > 0:
                print(f"  {bins[-1]:5d}+     : {count_large:6d} runs")

            # Show the 20 largest runs
            all_runs_sorted = sorted(all_runs, key=lambda x: -x[2])
            print(f"\nLargest interior invalid runs:")
            print(f"  {'Row':>5s}  {'ColStart':>8s}  {'Length':>7s}  {'ColEnd':>8s}")
            for r, cs, l in all_runs_sorted[:30]:
                print(f"  {r:5d}  {cs:8d}  {l:7d}  {cs+l-1:8d}")

        # ═══════ C: Check if interior invalids form vertical stripes ═══════
        print("\n" + "=" * 72)
        print("C. VERTICAL STRIPE ANALYSIS (interior columns)")
        print("=" * 72)
        # Check if certain columns are consistently invalid across many rows
        interior_inv_rate = inv_col_histogram / n_rows_sampled
        high_inv_cols = np.where(interior_inv_rate > 0.5)[0]
        if len(high_inv_cols) > 0:
            print(f"Columns with >50% invalid in rows 3000–7500: {len(high_inv_cols)}")
            # Find contiguous groups
            groups = []
            if len(high_inv_cols) > 0:
                grp_start = high_inv_cols[0]
                prev = high_inv_cols[0]
                for c in high_inv_cols[1:]:
                    if c > prev + 1:
                        groups.append((grp_start, prev - grp_start + 1))
                        grp_start = c
                    prev = c
                groups.append((grp_start, prev - grp_start + 1))

            print(f"Contiguous high-invalid column groups ({len(groups)}):")
            for gs, gl in groups[:20]:
                print(f"  cols {gs}–{gs+gl-1}  (width={gl})")
        else:
            print("No columns with >50% invalid in interior rows.")

        # Also check moderate invalid rate columns
        mod_inv_cols = np.where((interior_inv_rate > 0.2) & (interior_inv_rate <= 0.5))[0]
        if len(mod_inv_cols) > 0:
            print(f"\nColumns with 20–50% invalid in rows 3000–7500: {len(mod_inv_cols)}")
            groups2 = []
            if len(mod_inv_cols) > 0:
                grp_start = mod_inv_cols[0]
                prev = mod_inv_cols[0]
                for c in mod_inv_cols[1:]:
                    if c > prev + 1:
                        groups2.append((grp_start, prev - grp_start + 1))
                        grp_start = c
                    prev = c
                groups2.append((grp_start, prev - grp_start + 1))

            print(f"Contiguous moderate-invalid column groups ({len(groups2)}):")
            for gs, gl in groups2[:20]:
                rate = np.mean(interior_inv_rate[gs:gs+gl])
                print(f"  cols {gs}–{gs+gl-1}  (width={gl}, mean inv rate={rate*100:.1f}%)")

        # ═══════ D: Row-level invalid within footprint ═══════
        print("\n" + "=" * 72)
        print("D. PER-ROW INTERIOR INVALID COUNT (rows 3000–7500, within footprint)")
        print("=" * 72)
        # For well-covered rows, count invalid WITHIN the valid range
        print(f"  {'Row':>5s}  {'FootprintW':>10s}  {'InvInFP':>8s}  {'InvFrac':>8s}")
        for r in range(3000, 7500, 100):
            v_cols = np.where(valid[r, :])[0]
            if len(v_cols) < 100:
                continue
            first_v = v_cols[0]
            last_v = v_cols[-1]
            footprint_w = last_v - first_v + 1
            inv_in_fp = int(np.sum(invalid[r, first_v:last_v+1]))
            frac = inv_in_fp / footprint_w if footprint_w > 0 else 0
            if frac > 0.01:  # Only show rows with >1% interior invalid
                print(f"  {r:5d}  {footprint_w:10d}  {inv_in_fp:8d}  {frac*100:7.1f}%")

        # ═══════ E: Check benchmark metadata consistency ═══════
        print("\n" + "=" * 72)
        print("E. METADATA VS BENCHMARK CROSS-CHECK")
        print("=" * 72)
        print(f"Raster shape:           {H} × {W}")
        print(f"Benchmark expected:     10078 × 12263")
        print(f"Match:                  {'YES' if H == 10078 and W == 12263 else 'NO'}")
        print(f"Valid pixels:           {int(np.sum(valid)):,}")
        print(f"Benchmark expected:     91,113,091")
        print(f"Match:                  {'YES' if int(np.sum(valid)) == 91113091 else 'NO — diff=' + str(int(np.sum(valid)) - 91113091)}")
        print(f"Valid %:                {100*np.sum(valid)/(H*W):.5f}%")
        print(f"Benchmark expected:     73.72414%")

        # Pixel spacing from transform
        print(f"\nPixel spacing X:        {abs(ds.transform.a):.5f} m")
        print(f"Pixel spacing Y:        {abs(ds.transform.e):.5f} m")
        print(f"Benchmark range sp:     11.64781 m")
        print(f"Benchmark azimuth sp:   13.89642 m")
        print(f"Benchmark multilook:    range=5, azimuth=1")
        print(f"Note: 20m output != native multilook spacing — geocoding resampled to 20m grid")

    print("\nDone.")


if __name__ == "__main__":
    main()
