#!/bin/bash
# Quick test to verify subswath metadata integration

set -e

export RUST_LOG="info,sardine=debug"
export SARDINE_SERDE_ONLY=1
export SARDINE_REQUIRE_SUBSWATHS=1

cd /home/datacube/apps/SARdine

echo "==========================================================================="
echo "🧪 Quick Subswath Metadata Integration Test"
echo "==========================================================================="
echo ""

# Rebuild first
echo "📦 Rebuilding with latest changes..."
cd SARdine
pip install -e . --no-build-isolation 2>&1 | tail -10
cd ..

echo ""
echo "🚀 Running test..."
echo ""

python3 -m sardine.cli backscatter \
    data/S1A_IW_SLC__1SDV_20201230T165244_20201230T165311_035918_0434F0_6788.SAFE \
    pipeline_output/quick_subswath_test \
    --polarization VV \
    --optimization-mode complete \
    --range-looks 2 \
    --azimuth-looks 1 \
    --no-speckle-filter \
    --sequential \
    2>&1 | tee pipeline_output/quick_subswath_test.log | grep -E "Formatting.*subswath|Loaded.*subswath|Step 11|valid pixels|Range-Doppler|EPOCH|complete" | head -50

echo ""
echo "==========================================================================="
echo "📊 Results"
echo "==========================================================================="

# Check for subswath loading
echo ""
echo "Subswath loading:"
grep -E "Formatting.*subswath|Loaded.*subswath geometries" pipeline_output/quick_subswath_test.log | tail -5

# Check for geocoding results
echo ""
echo "Geocoding results:"
grep -E "valid pixels:|Geocoding|terrain correction" pipeline_output/quick_subswath_test.log | tail -10

echo ""
echo "Full log: pipeline_output/quick_subswath_test.log"
