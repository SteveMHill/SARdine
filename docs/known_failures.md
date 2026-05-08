# Known Failures

## Confirmed recurring symptoms
- Seams / stripes between bursts
- Burst duplication or overlap problems
- Unstable geocoding / terrain correction
- High coupling: fixing one thing often breaks another part of the pipeline

## Working assumptions
These symptoms may indicate a mix of:
- algorithmic defects
- indexing / overlap logic defects
- geometry / coordinate handling defects
- orchestration / state propagation defects

## Important rule
We do not assume the currently implemented pipeline stages are correct just because they exist.

## Case 001 observed artifact pattern
For input `S1A_IW_SLC__1SDV_20201005T170824_20201005T170851_034664_04098A_1E66.SAFE` in VV / sigma0 / dB output:
- final geocoded product shows duplicated-looking upper strips
- final product contains long blank band-like missing regions
- output metadata exists
- screenshots exist
- current interpretation: likely structural processing defect, not simple visualization issue