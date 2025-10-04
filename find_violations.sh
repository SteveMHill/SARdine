#!/bin/bash
# Find refactoring violations in SARdine codebase
# Usage: ./find_violations.sh

set -e

cd "$(dirname "$0")/SARdine"

echo "=================================================="
echo "🔍 SEARCHING FOR REFACTORING VIOLATIONS"
echo "=================================================="
echo ""

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

violations_found=0

# 1. Incidence angle hacks
echo -e "${YELLOW}1. Incidence Angle Hacks (should only be in IncidenceAngleModel)${NC}"
echo "   Files: src/io/slc_reader.rs, src/core/deburst.rs"
if rg -n "incidence_angle.*=.*(20\.0|linear)|min_incidence.*max_incidence|INCIDENCE_ANGLE_MIN|INCIDENCE_ANGLE_MAX" \
   src/io/slc_reader.rs src/core/deburst.rs 2>/dev/null; then
    echo -e "${RED}   ❌ VIOLATION FOUND${NC}"
    violations_found=$((violations_found + 1))
else
    echo -e "${GREEN}   ✅ Clean${NC}"
fi
echo ""

# 2. Calibration functions in wrong places
echo -e "${YELLOW}2. Calibration Functions (should only be in calibrate.rs)${NC}"
echo "   Files: src/io/*.rs, src/core/deburst.rs"
if rg -n "fn.*calibrate|calibrate_.*burst|apply.*calibration" \
   src/io/ src/core/deburst.rs 2>/dev/null | grep -v "parse_calibration_from_xml"; then
    echo -e "${RED}   ❌ VIOLATION FOUND${NC}"
    violations_found=$((violations_found + 1))
else
    echo -e "${GREEN}   ✅ Clean${NC}"
fi
echo ""

# 3. Noise removal in wrong places
echo -e "${YELLOW}3. Noise Removal (should only be in calibrate.rs)${NC}"
echo "   Files: src/io/*.rs, src/core/deburst.rs"
if rg -n "fn.*noise.*removal|subtract.*noise|apply.*noise.*removal" \
   src/io/ src/core/deburst.rs 2>/dev/null | grep -v "parse_noise_from_xml" | grep -v "find_noise_files" | grep -v "is_noise_xml"; then
    echo -e "${RED}   ❌ VIOLATION FOUND${NC}"
    violations_found=$((violations_found + 1))
else
    echo -e "${GREEN}   ✅ Clean${NC}"
fi
echo ""

# 4. Antenna pattern application in wrong places
echo -e "${YELLOW}4. Antenna Pattern Application (should only be in calibrate.rs)${NC}"
echo "   Files: src/io/*.rs, src/core/deburst.rs"
if rg -n "apply.*antenna|antenna.*correction" \
   src/io/ src/core/deburst.rs 2>/dev/null | grep -v "parse_antenna_pattern" | grep -v "get_antenna_patterns"; then
    echo -e "${RED}   ❌ VIOLATION FOUND${NC}"
    violations_found=$((violations_found + 1))
else
    echo -e "${GREEN}   ✅ Clean${NC}"
fi
echo ""

# 5. Dense LUT building in parsers
echo -e "${YELLOW}5. Dense LUT Building in Parsers (should only be in calibrate.rs)${NC}"
echo "   File: src/io/annotation.rs"
if rg -n "vec!\[.*[;,].*width\]|Array2::zeros.*\(.*width.*height|precompute.*lut" \
   src/io/annotation.rs 2>/dev/null | grep -v "// Example" | grep -v "width: usize"; then
    echo -e "${RED}   ❌ VIOLATION FOUND${NC}"
    violations_found=$((violations_found + 1))
else
    echo -e "${GREEN}   ✅ Clean${NC}"
fi
echo ""

# 6. Coordinate mapping in wrong places
echo -e "${YELLOW}6. Coordinate Mapping (should only be in calibrate.rs)${NC}"
echo "   Files: src/io/slc_reader.rs, src/core/deburst.rs"
if rg -n "pixel.*to.*lut|map.*coordinate|remap.*pixel|CalibrationCoordinateMapper" \
   src/io/slc_reader.rs src/core/deburst.rs 2>/dev/null | grep -v "// " | grep -v "///"; then
    echo -e "${RED}   ❌ VIOLATION FOUND${NC}"
    violations_found=$((violations_found + 1))
else
    echo -e "${GREEN}   ✅ Clean${NC}"
fi
echo ""

# 7. Power/intensity conversion in wrong places (SAR-specific, not geometric calculations)
echo -e "${YELLOW}7. Power/Intensity Conversion (should only be in calibrate.rs)${NC}"
echo "   Files: src/io/*.rs, src/core/deburst.rs"
# Exclude false positives: velocity_magnitude, gradient_magnitude (legitimate geometric calculations)
if rg -n "norm_sqr|magnitude.*sqr|to_power|complex.*intensity|real.*imag.*intensity" \
   src/io/ src/core/deburst.rs 2>/dev/null | \
   grep -v "// " | \
   grep -v "velocity_magnitude" | \
   grep -v "gradient_magnitude"; then
    echo -e "${RED}   ❌ VIOLATION FOUND${NC}"
    violations_found=$((violations_found + 1))
else
    echo -e "${GREEN}   ✅ Clean${NC}"
fi
echo ""

# Summary
echo "=================================================="
if [ $violations_found -eq 0 ]; then
    echo -e "${GREEN}✅ ALL CHECKS PASSED - NO VIOLATIONS FOUND${NC}"
    exit 0
else
    echo -e "${RED}❌ FOUND $violations_found VIOLATION(S)${NC}"
    echo ""
    echo "Next steps:"
    echo "  1. Review violations above"
    echo "  2. Create GitHub issues for each violation"
    echo "  3. Follow REFACTORING_CLEANUP_PLAN.md to fix"
    echo "  4. Re-run this script after fixes"
    exit 1
fi
