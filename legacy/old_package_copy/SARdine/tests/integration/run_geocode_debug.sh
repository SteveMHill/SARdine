#!/bin/bash
# Run geocoding test with full DEBUG logging to diagnose Range-Doppler failure

export RUST_LOG="debug,sardine::core::terrain_correction=debug"
export SARDINE_SERDE_ONLY=1
export SARDINE_REQUIRE_SUBSWATHS=1

cd /home/datacube/apps/SARdine

echo "==================================================================="
echo "🔍 Geocoding Debug Test with DEBUG Logging"
echo "==================================================================="
echo ""
echo "This will show WHY Range-Doppler transformation is failing."
echo "Looking for:"
echo "  - Newton-Raphson errors"
echo "  - Orbit interpolation errors"
echo "  - Epoch mismatch errors"
echo "  - Invalid coordinate errors"
echo ""

python3 -m sardine.cli backscatter \
    data/S1A_IW_SLC__1SDV_20201230T165244_20201230T165311_035918_0434F0_6788.SAFE \
    /tmp/geocode_debug_output \
    --polarization VV \
    --optimization-mode complete \
    --range-looks 2 \
    --azimuth-looks 1 \
    --no-speckle-filter \
    --sequential \
    2>&1 | tee /tmp/geocode_debug_full.log

echo ""
echo "==================================================================="
echo "📊 Analyzing failures..."
echo "==================================================================="

# Extract key error messages
echo ""
echo "Newton-Raphson failures:"
grep -E "Newton-Raphson failed|Newton-Raphson diverged|non-finite range_dot_velocity" /tmp/geocode_debug_full.log | head -10

echo ""
echo "Orbit interpolation failures:"
grep -E "Orbit interpolation failed" /tmp/geocode_debug_full.log | head -10

echo ""
echo "Epoch mismatch errors:"
grep -E "EPOCH MISMATCH" /tmp/geocode_debug_full.log | head -10

echo ""
echo "Invalid coordinate errors:"
grep -E "Invalid input coordinates|Invalid ECEF|Invalid slant range" /tmp/geocode_debug_full.log | head -10

echo ""
echo "Time windowing rejections:"
grep -E "Time windowing.*outside product window" /tmp/geocode_debug_full.log | head -10

echo ""
echo "==================================================================="
echo "Full debug log saved to: /tmp/geocode_debug_full.log"
echo "==================================================================="
