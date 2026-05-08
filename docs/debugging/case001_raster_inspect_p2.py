#!/usr/bin/env python3
"""
Case 001 deep artifact analysis — Part 2.
Focused on:
  A) Fine-grained invalid structure in the upper scene region (rows 0–2500)
  B) Horizontal stripe detection via row power variance
  C) Duplication detection via block cross-correlation
  D) Left/right edge invalid pattern (diagonal footprint shape)
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

        # ═══════ A: Fine-grained invalid profile in upper region ═══════
        print("=" * 72)
        print("A. FINE-GRAINED INVALID PROFILE — UPPER REGION (rows 0–2500)")
        print("=" * 72)
        # Per-row invalid fraction, every 10 rows
        print(f"  {'Row':>5s}  {'InvFrac':>7s}  {'ValidFirst':>10s}  {'ValidLast':>10s}  {'ValidSpan':>10s}")
        for r in range(0, min(2500, H), 25):
            row_inv = invalid[r, :]
            inv_frac = np.mean(row_inv)
            valid_cols = np.where(~row_inv)[0]
            if len(valid_cols) > 0:
                vf = valid_cols[0]
                vl = valid_cols[-1]
                vs = vl - vf + 1
            else:
                vf = vl = vs = -1
            print(f"  {r:5d}  {inv_frac*100:6.1f}%  {vf:10d}  {vl:10d}  {vs:10d}")

        # ═══════ B: Row power variance profile (stripe detection) ═══════
        print("\n" + "=" * 72)
        print("B. ROW-MEAN POWER PROFILE (fine-grained, rows 0–2500)")
        print("=" * 72)
        # For each row, compute mean of valid pixels
        row_means = np.full(H, np.nan)
        row_valid_counts = np.sum(valid, axis=1)
        for r in range(H):
            if row_valid_counts[r] > 100:
                row_means[r] = np.nanmean(data[r, valid[r, :]])

        print(f"  {'Row':>5s}  {'Mean(dB)':>9s}  {'ValidPx':>8s}  {'Diff':>7s}")
        prev_mean = np.nan
        for r in range(0, min(2500, H), 10):
            if np.isfinite(row_means[r]):
                diff = row_means[r] - prev_mean if np.isfinite(prev_mean) else 0.0
                print(f"  {r:5d}  {row_means[r]:9.3f}  {row_valid_counts[r]:8d}  {diff:+7.3f}")
                prev_mean = row_means[r]
            else:
                print(f"  {r:5d}       NaN  {row_valid_counts[r]:8d}")

        # ═══════ C: Row-mean discontinuities (jumps > 1 dB) ═══════
        print("\n" + "=" * 72)
        print("C. ROW-MEAN DISCONTINUITIES (jumps > 0.5 dB)")
        print("=" * 72)
        print(f"  {'Row':>5s}  {'Mean':>9s}  {'PrevMean':>9s}  {'Jump':>7s}")
        prev = np.nan
        jump_count = 0
        for r in range(H):
            if np.isfinite(row_means[r]) and np.isfinite(prev):
                jump = row_means[r] - prev
                if abs(jump) > 0.5:
                    print(f"  {r:5d}  {row_means[r]:9.3f}  {prev:9.3f}  {jump:+7.3f}")
                    jump_count += 1
            prev = row_means[r] if np.isfinite(row_means[r]) else prev
        print(f"\nTotal jumps > 0.5 dB: {jump_count}")

        # ═══════ D: Block duplication detection ═══════
        print("\n" + "=" * 72)
        print("D. BLOCK DUPLICATION DETECTION")
        print("=" * 72)
        # Compare blocks of ~100 rows against every other block
        # Focus on areas with >80% valid pixels
        block_size = 100
        # Pick a reference column slice (center of image)
        col_start = W // 4
        col_end = 3 * W // 4
        blocks = []
        for r_start in range(0, H - block_size, block_size):
            r_end = r_start + block_size
            block = data[r_start:r_end, col_start:col_end].copy()
            block_valid = valid[r_start:r_end, col_start:col_end]
            valid_frac = np.mean(block_valid)
            if valid_frac > 0.5:
                # Set invalid to 0 for comparison purposes
                block[~block_valid] = 0
                blocks.append((r_start, block, block_valid, valid_frac))

        print(f"Comparing {len(blocks)} blocks of {block_size} rows "
              f"(cols {col_start}–{col_end})")

        # Find high-correlation block pairs (potential duplicates)
        high_corr_pairs = []
        for i in range(len(blocks)):
            for j in range(i + 1, min(i + 40, len(blocks))):  # Look within 4000 rows
                r_i, b_i, v_i, _ = blocks[i]
                r_j, b_j, v_j, _ = blocks[j]
                shared = v_i & v_j
                n_shared = np.sum(shared)
                if n_shared < block_size * (col_end - col_start) * 0.3:
                    continue
                a = b_i[shared].astype(np.float64)
                b = b_j[shared].astype(np.float64)
                a_c = a - np.mean(a)
                b_c = b - np.mean(b)
                denom = np.sqrt(np.sum(a_c**2) * np.sum(b_c**2))
                if denom == 0:
                    continue
                corr = np.sum(a_c * b_c) / denom
                if corr > 0.7:
                    mae = float(np.mean(np.abs(a - b)))
                    high_corr_pairs.append((r_i, r_j, corr, mae, n_shared))

        if high_corr_pairs:
            # Sort by correlation descending
            high_corr_pairs.sort(key=lambda x: -x[2])
            print(f"\nHigh-correlation block pairs (corr > 0.7): {len(high_corr_pairs)}")
            print(f"  {'Block1':>7s}  {'Block2':>7s}  {'Offset':>7s}  {'Corr':>7s}  {'MAE':>9s}  {'SharedPx':>9s}")
            for r_i, r_j, corr, mae, nsh in high_corr_pairs[:30]:
                offset = r_j - r_i
                print(f"  {r_i:7d}  {r_j:7d}  {offset:7d}  {corr:7.4f}  {mae:9.4f}  {nsh:9d}")

            # Check if there's a dominant offset among duplicates
            offsets = [r_j - r_i for r_i, r_j, c, m, n in high_corr_pairs]
            unique_offsets, counts = np.unique(offsets, return_counts=True)
            print("\nOffset frequency among high-corr pairs:")
            for off, cnt in sorted(zip(counts, unique_offsets), reverse=True)[:10]:
                print(f"  offset={cnt:4d} rows:  {off} occurrences")
        else:
            print("\nNo high-correlation block pairs found (no obvious duplication).")

        # ═══════ E: Diagonal footprint shape analysis ═══════
        print("\n" + "=" * 72)
        print("E. VALID PIXEL FOOTPRINT SHAPE")
        print("=" * 72)
        # For each row, find first and last valid column
        print(f"  {'Row':>5s}  {'FirstCol':>8s}  {'LastCol':>8s}  {'Width':>6s}")
        for r in range(0, H, max(1, H // 30)):
            v_cols = np.where(valid[r, :])[0]
            if len(v_cols) > 0:
                print(f"  {r:5d}  {v_cols[0]:8d}  {v_cols[-1]:8d}  {v_cols[-1]-v_cols[0]+1:6d}")
            else:
                print(f"  {r:5d}      none      none       0")

    print("\nDone.")


if __name__ == "__main__":
    main()
