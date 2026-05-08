# Rebuild Boundary 01

## Purpose
Define what the rebuild is trying to achieve, and what legacy debugging is no longer meant to do.

## The rebuild is for
- building a new verified processing core
- making assumptions explicit
- enforcing contracts and invariants
- reducing hidden state and unsafe coupling
- enabling stage-by-stage verification
- selectively reusing only trustworthy legacy components

## The rebuild is not for
- fixing every bug in the legacy package
- fully understanding every legacy code path
- making the existing architecture reliable
- preserving old module boundaries just because they already exist

## Legacy code will be used for
- reference
- trust classification
- selective reuse
- comparison artifacts
- extraction of hidden assumptions

## Legacy code will not be used for
- open-ended debugging
- broad in-place repair
- gradual patching of the whole old pipeline

## Current conclusion
We have enough evidence to stop broad legacy debugging and move into rebuild design.

## Rebuild principle
The new pipeline becomes the source of truth.
The legacy pipeline remains a reference artifact.