"""Compare sigma_nought vs gamma LUT values in calibration XMLs."""
import xml.etree.ElementTree as ET
import glob
import math

SAFE = (
    "data/SLC/S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE"
)

for iw in ["iw1", "iw2", "iw3"]:
    cal_path = sorted(
        glob.glob(f"{SAFE}/annotation/calibration/calibration-s1b-{iw}-slc-vv-*.xml")
    )[0]
    root = ET.parse(cal_path).getroot()
    vecs = root.findall(".//calibrationVector")
    v0 = vecs[len(vecs)//2]  # middle vector

    px = list(map(int, v0.find("pixel").text.split()))
    sigma = list(map(float, v0.find("sigmaNought").text.split()))
    gamma = list(map(float, v0.find("gamma").text.split()))
    beta = list(map(float, v0.find("betaNought").text.split()))

    # Ratio at several sample points
    mid = len(px) // 2
    for idx in [0, mid, -1]:
        s = sigma[idx]
        g = gamma[idx]
        b = beta[idx]
        ratio_gs = g / s
        ratio_sb = s / b
        # If gamma = sigma / sqrt(cos(theta)), then gamma/sigma = 1/sqrt(cos(theta))
        # So cos(theta) = (sigma/gamma)^2
        cos_theta = (s / g) ** 2
        theta_deg = math.degrees(math.acos(cos_theta)) if -1 <= cos_theta <= 1 else float("nan")
        print(
            f"{iw.upper()} pix={px[idx]:5d}: "
            f"sigma={s:.3f}  gamma={g:.3f}  beta={b:.3f}  "
            f"gamma/sigma={ratio_gs:.5f}  implied_theta={theta_deg:.2f}°"
        )
    print()
