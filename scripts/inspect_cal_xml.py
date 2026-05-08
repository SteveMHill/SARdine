"""Inspect calibration LUT pixel coords and burst valid sample ranges."""
import xml.etree.ElementTree as ET
import glob

SAFE = (
    "data/SLC/S1B_IW_SLC__1SDV_20190123T053348_20190123T053415_014617_01B3D4_E833.SAFE"
)

for iw in ["iw1", "iw2", "iw3"]:
    cal_path = sorted(
        glob.glob(f"{SAFE}/annotation/calibration/calibration-s1b-{iw}-slc-vv-*.xml")
    )[0]
    root = ET.parse(cal_path).getroot()
    vecs = root.findall(".//calibrationVector")
    v0 = vecs[0]
    px = list(map(int, v0.find("pixel").text.split()))
    sigma = list(map(float, v0.find("sigmaNought").text.split()))
    print(f"{iw.upper()}: {len(vecs)} vectors, pixels[0]={px[0]}, pixels[-1]={px[-1]}, n_pix={len(px)}")
    print(f"  sigma_nought[0]={sigma[0]:.3f}  mid={sigma[len(sigma)//2]:.3f}  last={sigma[-1]:.3f}")

    ann_path = sorted(glob.glob(f"{SAFE}/annotation/s1b-{iw}-slc-vv-*.xml"))[0]
    aroot = ET.parse(ann_path).getroot()
    b0 = aroot.findall(".//burst")[0]
    fvs = list(map(int, b0.find("firstValidSample").text.strip().split()))
    lvs = list(map(int, b0.find("lastValidSample").text.strip().split()))
    print(f"  burst0 firstValidSample: min={min(fvs)}, max={max(fvs)}")
    print(f"  burst0 lastValidSample:  min={min(lvs)}, max={max(lvs)}")

    img_info = aroot.find(".//imageAnnotation/imageInformation")
    n_samples = int(img_info.find("numberOfSamples").text)
    n_lines = int(img_info.find("numberOfLines").text)
    print(f"  TIFF dims: {n_lines} lines x {n_samples} samples")
    print()
