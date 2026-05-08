# Rebuild Plan

## Goal
Rescue the SARdine codebase through controlled re-architecture, not blind rewrite.

## Phase order
1. Repository forensics
2. Failure diagnosis
3. Target architecture
4. Module contracts
5. Verification strategy
6. Minimal vertical slice
7. Controlled module-by-module replacement

## Rules
- Narrow scope aggressively
- Prefer explicit uncertainty over guessed implementations
- No placeholder production logic
- No broad refactors without justification
- Every important change must be verifiable