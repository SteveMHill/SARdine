# Legacy Status

## Purpose
This repository contains an existing Sentinel-1 SAR processing implementation with valuable domain knowledge and partial working code.

## Important
The current implementation is NOT treated as ground truth.

## Default rule
Legacy code must not be broadly refactored or rewritten in place.

## Allowed use of legacy code
Legacy code may be:
- inspected
- documented
- wrapped
- reused selectively if isolated and justified
- compared against new implementations

## Not allowed by default
- broad cleanup refactors
- style-driven rewrites
- replacing large sections without a contract
- treating existing behavior as correct without verification

## Required classification for reused legacy components
Every reused component must be explicitly classified as:
- safe to reuse
- reusable with modification
- unsafe / reference only

## Migration intent
We will gradually replace or wrap legacy functionality through contract-driven modules and verification, starting with the highest-value diagnostic areas.