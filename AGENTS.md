# AI working rules for this repository

This repository contains a legacy Sentinel-1 SAR backscatter processing pipeline in Python and Rust.

The legacy codebase contains useful domain knowledge, parsing logic, and implementation work.
However, parts of the pipeline are known to be incorrect or unstable.
Therefore:

## Core rule
Legacy code is reference material, not ground truth.

## Non-negotiable behavior rules
- Do not rewrite the whole pipeline.
- Do not make broad refactors unless explicitly requested.
- Do not introduce placeholder implementations, fake adapters, mocked success paths, or silent fallback behavior.
- Do not invent domain logic just to make the workflow appear complete.
- If correctness is uncertain, stop and state the uncertainty explicitly.
- Prefer a small explicit failure over a large plausible-but-wrong implementation.
- Reuse legacy code only if it is isolated, understandable, testable, and justified.
- Every reused legacy component must be classified as:
  - safe to reuse
  - reusable with modification
  - unsafe / reference only

## Working mode
For non-trivial tasks:
1. Inspect relevant files first.
2. Summarize findings.
3. List assumptions.
4. List risks and uncertainties.
5. Propose the smallest safe next step.
6. Only then edit code.

## Output requirements
For every meaningful task, provide:
1. What you inspected
2. What you think is true
3. What you are unsure about
4. Smallest safe next change
5. How it will be verified

## Verification rules
No logic change is complete until:
- relevant tests are added or updated
- the narrowest relevant test set is run
- verification results are summarized
- remaining uncertainty is stated explicitly

## Domain-specific SAR rules
Be explicit about:
- timing domains
- slant-range vs ground-range assumptions
- burst overlap handling
- deburst assumptions
- radiometric domain transitions
- geocoding target grid definition
- nodata and invalid mask propagation

## Forbidden shortcuts
- No “TODO-backed” production logic
- No fake completeness
- No compile-only success claims
- No hidden scope expansion
- No changing unrelated files to make the task easier

## Mechanical guard against silent fallbacks and hardcoded magic numbers

A class of bugs has repeatedly corrupted the scientific value of this package:
silent fallbacks (`unwrap_or`, `unwrap_or_else`, `unwrap_or_default`, `.ok()?`,
`let _ = result`) that swallow physically meaningful failures, and hardcoded
Sentinel-1 constants (azimuth time interval, pixel spacing, burst cycle,
carrier frequency, lines-per-burst) that should always come from annotation
metadata.

To prevent regressions, `scripts/check_no_silent_fallbacks.sh` is enforced as
an integration test (`new_pipeline/scene/tests/no_silent_fallbacks.rs`) and
runs on every `cargo test`.  The script blocks:

| Pattern class           | Examples                                          |
|-------------------------|---------------------------------------------------|
| Silent fallbacks        | `.unwrap_or(`, `.unwrap_or_else(`, `.unwrap_or_default(`, `.ok()?`, `let _ = …` |
| Unfinished logic        | `todo!`, `unimplemented!`, `panic!`, `// TODO/FIXME/XXX` |
| Hardcoded S-1 constants | ATI ≈ 2.05 ms, az pixel spacing ≈ 13–14 m, burst cycle ≈ 2.75–2.76 s, carrier 5.405 GHz, lines-per-burst ≈ 14xx |

### Two ways to satisfy the guard

1. **Replace** the pattern with an explicit error path on a typed error enum,
   or plumb the value from metadata.  This is the default expectation.
2. **Annotate** with a same-line `// SAFETY-OK: <reason>` trailing comment,
   stating the precise reason the use is provably safe (e.g. a chrono
   microseconds extraction over an orbit-window duration cannot overflow,
   or a `strip_prefix` whose no-match return value is the design intent).

The annotation must be on the *same physical line* as the offending token.
Block comments and prior-line comments do not count.  This forces the
intent into code review.

### PR self-checklist (mandatory)

Copy this into every PR description that touches `new_pipeline/scene/src/`:

- [ ] No new `unwrap_or*`, `unwrap_or_default`, `let _ = result`, `.ok()?` patterns (or each is annotated `SAFETY-OK: …` on the same line)
- [ ] No new hardcoded Sentinel-1 constants in production code (read from `SubSwathMetadata` / annotation parser instead)
- [ ] Every new failure mode is an enum variant on a typed error, not a logged warning or stderr print
- [ ] Every new "fallback" path is gated behind an explicit env var or function arg the caller must set
- [ ] At least one negative test asserts the new error variant fires
- [ ] `cargo test` passes (which runs `no_silent_fallbacks` integration test)

### Local verification

```sh
./scripts/check_no_silent_fallbacks.sh    # exits 0 if clean, 1 if violations
cargo test --test no_silent_fallbacks     # same check via cargo
cargo test                                # full suite, includes the guard
```