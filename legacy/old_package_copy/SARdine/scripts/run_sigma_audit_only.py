#!/usr/bin/env python3
"""Run Sigma0 calibration only and trigger sigma_audit JSON writing.

This script avoids GeoTIFF export and focuses on producing sigma_audit_*.json
for each IW subswath of the chosen polarization.
"""

import os
import sys
import time
import logging
from pathlib import Path
from argparse import ArgumentParser

# Set logging and Rust env before importing sardine
os.environ.setdefault("RUST_LOG", "info")
os.environ.setdefault("PYTHONUNBUFFERED", "1")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


def parse_args() -> ArgumentParser:
    p = ArgumentParser(description="Run Sigma0 calibration and emit sigma_audit JSONs only")
    p.add_argument("--product", required=True, help="Path to Sentinel-1 SAFE product")
    p.add_argument("--polarization", "-p", default="VV", help="Polarization (e.g. VV, VH)")
    p.add_argument(
        "--subswath",
        "-s",
        default="ALL",
        help="Specific IW subswath (e.g. IW2) or ALL for all available",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()

    safe_path = Path(args.product).expanduser().resolve()
    pol = args.polarization.upper()
    subswath_filter = args.subswath.upper()

    if not safe_path.exists():
        logger.error(f"SAFE product not found: {safe_path}")
        return 1

    audit_dir = os.environ.get("SARDINE_CALIB_AUDIT_DIR", "")
    if not audit_dir:
        logger.warning(
            "SARDINE_CALIB_AUDIT_DIR not set – no sigma audit JSONs will be written."
        )
    else:
        logger.info(f"Sigma audit directory: {audit_dir}")

    import sardine

    logger.info(f"Opening SAFE: {safe_path}")
    reader = sardine.create_cached_slc_reader(str(safe_path))
    metadata = reader.get_cached_metadata()

    subswaths_by_pol = reader.get_all_iw_subswaths()
    if pol not in subswaths_by_pol:
        logger.error(f"Polarization {pol} not found in product.")
        return 1

    subswaths_dict = subswaths_by_pol[pol]
    all_subswaths = sorted(subswaths_dict.keys())

    if subswath_filter != "ALL":
        if subswath_filter not in subswaths_dict:
            logger.error(
                f"Requested subswath {subswath_filter} not found for {pol}. "
                f"Available: {all_subswaths}"
            )
            return 1
        target_subswaths = [subswath_filter]
    else:
        target_subswaths = all_subswaths

    logger.info(
        f"Will run Sigma0 calibration (no exports) for {pol} on subswaths: {target_subswaths}"
    )

    total_start = time.time()

    for sw_name in target_subswaths:
        logger.info("=" * 60)
        logger.info(f"Sigma0 calibration for {sw_name} {pol}")
        step_start = time.time()

        try:
            deburst_result = sardine.deburst_topsar_cached(reader, sw_name, pol)

            # unwrap dict-style result if present
            if isinstance(deburst_result, dict):
                if deburst_result.get("status") == "error":
                    logger.error(f"  Deburst error for {sw_name}: {deburst_result.get('message')}")
                    continue
                if "data" in deburst_result:
                    slc_data = deburst_result["data"]
                else:
                    slc_data = None
                    for k, v in deburst_result.items():
                        import numpy as np

                        if isinstance(v, np.ndarray) and v.ndim >= 2:
                            slc_data = v
                            logger.info(f"  Using deburst data from key '{k}'")
                            break
                    if slc_data is None:
                        logger.error("  No ndarray payload found in deburst result; skipping.")
                        continue
            else:
                slc_data = deburst_result

            logger.info(f"  Deburst data shape: {getattr(slc_data, 'shape', '?')} ")

            # Use cached-reader calibration with thermal noise removal enabled.
            # This path also triggers sigma-audit JSON writing when
            # SARDINE_CALIB_AUDIT_DIR is set in the environment.
            calib_result = sardine.radiometric_calibration_with_denoising_cached(
                reader,
                sw_name,
                pol,
                "sigma0",
                slc_data,
                True,
            )

            if isinstance(calib_result, dict) and calib_result.get("status") == "error":
                logger.error(f"  Calibration error for {sw_name}: {calib_result.get('message')}")
                continue

            logger.info(
                f"  Sigma0 calibration for {sw_name} completed in "
                f"{time.time() - step_start:.2f}s (array discarded; audit should be written)."
            )

        except Exception as e:  # noqa: BLE001
            logger.error(f"  Exception during calibration for {sw_name}: {e}")
            import traceback

            traceback.print_exc()

    logger.info("=" * 60)
    logger.info(f"Total elapsed: {time.time() - total_start:.2f}s")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
