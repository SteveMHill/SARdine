# Contributing to SARdine

Thanks for your interest in improving SARdine. This document captures the
non-negotiable rules and the local workflow.

## Read these first

1. [`AGENTS.md`](AGENTS.md) — the working rules. Legacy code is reference
   material, not ground truth. No silent fallbacks. No hardcoded
   Sentinel-1 constants.
2. [`LEGACY_STATUS.md`](LEGACY_STATUS.md) — what is allowed and forbidden
   when touching `legacy/`.
3. [`docs/architecture_overview.md`](docs/architecture_overview.md) — the
   current pipeline shape.

## Local development

```sh
# build + test (runs the no-silent-fallbacks guard automatically)
cargo test --workspace

# guard alone
./scripts/check_no_silent_fallbacks.sh

# formatting and lint
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
```

The Python extension is built via `maturin develop` from
`sardine-py/`.

## Pull-request checklist

Copy this into every PR description that touches `sardine/src/`:

- [ ] No new `unwrap_or*`, `unwrap_or_default`, `let _ = result`, `.ok()?`
      patterns (or each is annotated `// SAFETY-OK: …` on the same line).
- [ ] No new hardcoded Sentinel-1 constants in production code (read from
      `SubSwathMetadata` / annotation parser instead).
- [ ] Every new failure mode is an enum variant on a typed error, not a
      logged warning or stderr print.
- [ ] Every new "fallback" path is gated behind an explicit env var or
      function arg the caller must set.
- [ ] At least one negative test asserts the new error variant fires.
- [ ] `cargo test` passes (which runs the `no_silent_fallbacks`
      integration test).
- [ ] `CHANGELOG.md` updated under the `## [Unreleased]` section if the
      change is user-visible.

## Reporting bugs

Please include:
- The exact `sardine` invocation (CLI flags or Python call).
- Inputs: SAFE product ID, orbit file (POEORB/RESORB), DEM source, geoid.
- Output: stderr log (`sardine` writes timings and counters to stderr),
  the `.provenance.json` sidecar, and a small region of the failing raster
  if relevant.
- Expected vs observed behaviour.

## Code review

Numerical correctness is the highest bar. A change that ships faster code
at the cost of unverified output will be rejected. Phase-3 of the perf
work is a worked example of what "bit-equivalent on a real scene" looks
like — see the regression tests in `sardine/src/orbit.rs`.
