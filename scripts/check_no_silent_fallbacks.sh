#!/usr/bin/env bash
# Forbidden-pattern scanner for the new_pipeline/scene crate.
#
# We have been bitten repeatedly by "silent fallback" patterns that swallow
# physically meaningful failures and replace them with plausible-but-wrong
# defaults.  This script enforces an opt-in policy: any use of a forbidden
# pattern must be annotated on the same line with a trailing comment of the
# form `// SAFETY-OK: <reason>` to be accepted.
#
# Exits 0 if clean, 1 if any unannotated offender is found.
#
# Run from repo root:  ./scripts/check_no_silent_fallbacks.sh

set -u

SRC_DIR="${1:-new_pipeline/scene/src}"

if [ ! -d "$SRC_DIR" ]; then
    echo "error: source directory not found: $SRC_DIR" >&2
    exit 2
fi

# Patterns we care about.  Each entry: "label|kind|extended-regex".
#   kind = strict   → match everywhere in src/
#   kind = nondoc   → skip doc comments, assert*! lines, and lines whose
#                     trimmed prefix is `//` (regular comments).  Used for
#                     magic-number patterns where reference values inside
#                     tests are legitimate but production-code literals are
#                     a silent-fallback hazard.
PATTERNS=(
    "unwrap_or         |strict|\\.unwrap_or\\("
    "unwrap_or_else    |strict|\\.unwrap_or_else\\("
    "unwrap_or_default |strict|\\.unwrap_or_default\\("
    "ok_questionmark   |strict|\\.ok\\(\\)\\?"
    "let_underscore    |strict|^[[:space:]]*let _[[:space:]]*=[[:space:]]"
    "todo_macro        |strict|^[[:space:]]*todo!\\("
    "unimplemented     |strict|^[[:space:]]*unimplemented!\\("
    "panic_macro       |strict|^[[:space:]]*panic!\\("
    "todo_comment      |strict|//[[:space:]]*(TODO|FIXME|XXX)"
    # ── Banned Sentinel-1 magic numbers (production code only) ─────────
    # These constants vary per scene and MUST come from annotation metadata.
    # Hardcoding them is the silent-fallback class of bug we want to prevent.
    # Reference values inside test assertions and doc comments are allowed.
    "magic_ati         |nondoc|0\\.00205[0-9]+|2\\.05[0-9]+e-?0?3"
    "magic_az_spacing  |nondoc|\\b1[34]\\.[0-9]+(_?f6?4)?\\b"
    "magic_burst_cycle |nondoc|2\\.7[5-6][0-9]+(_?f6?4)?"
    "magic_wavelength  |nondoc|0\\.0555[0-9]+|0\\.0556[0-9]*"
    "magic_carrier     |nondoc|5\\.405e\\+?0?9|5405000000\\.0"
    "magic_lpb         |nondoc|\\b14[0-9][0-9]_?usize\\b"
)

# Exit-zero accumulator
status=0

# We deliberately walk only .rs files under src/ (tests/, examples/, benches/
# are out of scope; #[cfg(test)] blocks inside src/ are still scanned because
# inline tests must hold the same standard for fixtures touching real data).
mapfile -t FILES < <(find "$SRC_DIR" -type f -name '*.rs' | sort)

for f in "${FILES[@]}"; do
    # Detect when we are inside a `#[cfg(test)]` annotated module by tracking
    # brace depth.  We don't try to be a real Rust parser; we only need to
    # skip lines whose enclosing module is gated by cfg(test).  This is a
    # heuristic: it handles the common pattern `#[cfg(test)] mod tests { ... }`
    # at end-of-file but does NOT handle nested non-test modules within a
    # test module (acceptable trade-off; flag with SAFETY-OK if needed).
    in_test_mod_state=$(awk '
        BEGIN { in_test = 0; depth = 0; pending = 0 }
        /^[[:space:]]*#\[cfg\(test\)\]/ { pending = 1 }
        {
            line = $0
            # Per-line brace counting (string-blind; good enough for src/).
            opens  = gsub(/\{/, "&", line)
            closes = gsub(/\}/, "&", line)
            # Re-read $0 since gsub modified line var only via copy.
            line = $0
            opens  = gsub(/\{/, "&", line); line = $0
            closes = gsub(/\}/, "&", line); line = $0
            print (in_test ? "T" : "N")
            if (pending && opens > 0) { in_test = 1; depth = 0; pending = 0 }
            if (in_test) {
                depth += opens - closes
                if (depth <= 0 && opens == 0) { in_test = 0 }
            }
        }
    ' "$f")

    for entry in "${PATTERNS[@]}"; do
        IFS='|' read -r label kind regex <<< "$entry"
        label="$(echo "$label" | xargs)"
        kind="$(echo "$kind" | xargs)"

        while IFS= read -r hit; do
            line_no="${hit%%:*}"
            line_text="${hit#*:}"

            # Allow opt-in annotation on the same line.
            if [[ "$line_text" == *"SAFETY-OK:"* ]]; then
                continue
            fi

            # nondoc patterns: ignore doc comments, line comments, asserts,
            # and lines inside #[cfg(test)] modules.
            if [[ "$kind" == "nondoc" ]]; then
                trimmed="${line_text#"${line_text%%[![:space:]]*}"}"
                if [[ "$trimmed" == //* ]]; then
                    continue
                fi
                if [[ "$line_text" == *"assert!"* ]] \
                    || [[ "$line_text" == *"assert_eq!"* ]] \
                    || [[ "$line_text" == *"assert_ne!"* ]] \
                    || [[ "$line_text" == *"debug_assert"* ]]; then
                    continue
                fi
                # Lookup test-mod state for this line number.
                state=$(echo "$in_test_mod_state" | sed -n "${line_no}p")
                if [[ "$state" == "T" ]]; then
                    continue
                fi
            fi

            echo "FORBIDDEN[$label] $f:$line_no: $line_text"
            status=1
        done < <(grep -E -n "$regex" "$f" 2>/dev/null || true)
    done
done

if [ "$status" -ne 0 ]; then
    cat >&2 <<'EOF'

================================================================================
Forbidden patterns detected.

Each match above is a likely silent-fallback or unfinished-logic site.  Either:

  1. Replace it with an explicit error path on a typed error enum, or
  2. If the use is provably safe (e.g. chrono microseconds overflow at >290k
     years, or a string-parser default that does not affect numerics), append
     a same-line comment:

         foo.unwrap_or(0); // SAFETY-OK: chrono microseconds cannot overflow

The annotation makes the intent reviewable.  Bare offenders break the build.
================================================================================
EOF
fi

exit "$status"
