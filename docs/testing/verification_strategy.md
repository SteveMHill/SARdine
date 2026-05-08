# Verification Strategy

## Principle
No important logic change is accepted without verification.

## Required verification layers
- Unit tests
- Contract tests
- Differential tests where useful
- Intermediate diagnostics
- Visual debugging artifacts where useful

## Domain-specific focus
For SAR processing, verification should eventually include:
- burst boundary diagnostics
- overlap / duplication diagnostics
- geocoding diagnostics
- invalid mask propagation checks
- intermediate product inspection

## Rule
A pipeline stage is not trusted because it runs.
It is trusted only when assumptions, invariants, and outputs are verified.