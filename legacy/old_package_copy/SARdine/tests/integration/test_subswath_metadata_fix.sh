#!/bin/bash
# Test script for subswath metadata fix in geocoding

set -e

echo "==========================================================================="
echo "🧪 Testing Subswath Metadata Fix for Range-Doppler Geocoding"
echo "==========================================================================="
echo ""

# Set environment
export RUST_LOG="info,sardine::core::terrain_correction=debug"
export SARDINE_SERDE_ONLY=1
export SARDINE_REQUIRE_SUBSWATHS=1

cd /home/datacube/apps/SARdine

# Rebuild the Rust extension with the new lib.rs changes
echo "📦 Rebuilding Rust extension..."
cd SARdine
pip install -e . --no-build-isolation 2>&1 | tail -20
cd ..

echo ""
echo "✅ Rebuild complete"
echo ""
echo "🚀 Running full pipeline test with geocoding..."
echo "   Looking for:"
echo "   - ✅ Loaded N subswath geometries from Python metadata"
echo "   - 📡 Loaded subswath IW1/IW2/IW3 debug messages"
echo "   - Range-Doppler success rate > 50%"
echo ""

# Run the CLI pipeline with geocoding
python3 -m sardine.cli backscatter \
    data/S1A_IW_SLC__1SDV_20201230T165244_20201230T165311_035918_0434F0_6788.SAFE \
    pipeline_output/subswath_fix_test \
    --polarization VV \
    --optimization-mode complete \
    --range-looks 2 \
    --azimuth-looks 1 \
    --no-speckle-filter \
    --sequential \
    2>&1 | tee pipeline_output/subswath_fix_test.log

echo ""
echo "==========================================================================="
echo "📊 Analysis of Results"
echo "==========================================================================="

# Check if subswaths were loaded
SUBSWATH_COUNT=$(grep -c "Loaded subswath" pipeline_output/subswath_fix_test.log || echo "0")
echo ""
echo "Subswaths loaded: $SUBSWATH_COUNT"
if [ "$SUBSWATH_COUNT" -ge 3 ]; then
    echo "   ✅ All 3 IW subswaths loaded successfully"
else
    echo "   ❌ Expected 3 subswaths, got $SUBSWATH_COUNT"
fi

# Check geocoding success rate
if grep -q "valid pixels:" pipeline_output/subswath_fix_test.log; then
    VALID_PERCENT=$(grep "valid pixels:" pipeline_output/subswath_fix_test.log | tail -1 | grep -oP '\d+\.\d+%')
    echo ""
    echo "Geocoding valid pixels: $VALID_PERCENT"
    PERCENT_VALUE=$(echo $VALID_PERCENT | grep -oP '[\d.]+')
    if (( $(echo "$PERCENT_VALUE > 50.0" | bc -l) )); then
        echo "   ✅ Geocoding success rate improved!"
    elif (( $(echo "$PERCENT_VALUE > 10.0" | bc -l) )); then
        echo "   ⚠️  Some improvement, but still low coverage"
    else
        echo "   ❌ No improvement in geocoding"
    fi
else
    echo ""
    echo "   ⚠️  Could not find geocoding statistics"
fi

# Check for Range-Doppler failures
RD_FAILURES=$(grep -c "Range-Doppler transformation failed" pipeline_output/subswath_fix_test.log || echo "0")
echo ""
echo "Range-Doppler failures: $RD_FAILURES"
if [ "$RD_FAILURES" -lt 1000 ]; then
    echo "   ✅ Significantly fewer Range-Doppler failures"
elif [ "$RD_FAILURES" -lt 100000 ]; then
    echo "   ⚠️  Still many Range-Doppler failures"
else
    echo "   ❌ Massive Range-Doppler failure count"
fi

# List output files
echo ""
echo "Output files created:"
ls -lh pipeline_output/subswath_fix_test/*.tif 2>/dev/null || echo "   ⚠️  No GeoTIFF files created"

echo ""
echo "==========================================================================="
echo "Full log: pipeline_output/subswath_fix_test.log"
echo "==========================================================================="
