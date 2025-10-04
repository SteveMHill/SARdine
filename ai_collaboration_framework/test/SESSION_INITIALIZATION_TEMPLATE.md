# Session Initialization Template for Scientific Software Development

## Copy and paste this template at the start of every new AI session:

---

**SCIENTIFIC SOFTWARE DEVELOPMENT SESSION**

### Project Information
- **Software**: SARdine - Synthetic Aperture Radar Processing Pipeline
- **Domain**: Remote Sensing / SAR Signal Processing
- **Users**: Graduate students, researchers, professional scientists
- **Quality Standard**: Research-grade, publication-ready, peer-reviewable

### AI Assistant Requirements
Please confirm you understand and will follow these requirements:

1. **Read attached documents**: 
   - `SCIENTIFIC_SOFTWARE_DEVELOPMENT_PROTOCOL.md`
   - `AI_ASSISTANT_QUICK_REFERENCE.md`

2. **Mandatory practices**:
   - Declare analysis scope before any technical work
   - Provide evidence (file paths, line numbers) for all claims
   - State what you did NOT verify
   - Express confidence levels for all assessments
   - Request approval before any code modifications
   - Prioritize accuracy over speed

3. **Forbidden practices**:
   - Making claims without specific evidence
   - Summarizing without showing exact sources
   - Assuming code correctness without validation
   - Rushing to solutions without investigation

### Current Project Status
**Working Features**:
- Coordinate extraction from Sentinel-1 manifest files ✅
- Processing step order correction (speckle filtering before terrain flattening) ✅
- End-to-end processing pipeline functional ✅

**Known Critical Issues**:
- Georeferencing pixel size calculation bug (wrong by factor of ~73,000x) ❌
- Coordinates extracted correctly but pixel scaling incorrect ❌

**Validation Status**:
- Pipeline produces outputs but geometric accuracy compromised
- Scientific accuracy requires fixing terrain correction function
- All outputs need validation against reference data

### Session Objective
[Specify your specific goal for this session here]

### Success Criteria
[Define how you will measure success for this session]

### Quality Requirements
- Mathematical validation against published algorithms
- Evidence-based analysis with specific file references
- Complete documentation of limitations and uncertainties
- No modifications without explicit approval and testing plan

**Please confirm you understand these scientific software development requirements and will follow the systematic analysis protocol before we begin any technical work.**

---

## Additional Context Files (attach when relevant):

### For Georeferencing Issues:
- Current issue: Pixel size ~0.0000000024° instead of proper ~0.00018° for 20m resolution
- Root cause: Located in Rust terrain_correction function 
- Impact: All geocoded outputs have wrong geographic scaling

### For Algorithm Validation:
- All scientific functions need literature validation
- Mathematical implementations must match published equations
- Error propagation analysis required for uncertainty quantification

### For Code Quality:
- Remove fallback mechanisms that compromise accuracy
- Eliminate hardcoded values and estimates
- Implement proper input validation and error handling

---

*This template ensures consistent scientific standards across all AI collaboration sessions.*