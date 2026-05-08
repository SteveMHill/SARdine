"""Strict science and validation wiring for BackscatterProcessor.

This module centralizes how "strict science" behaviour is derived from
options and environment variables and how it is propagated into the
process-wide environment. Keeping this logic here makes the main
processor implementation smaller and easier to read.
"""

from __future__ import annotations

from typing import Any, Dict, Tuple
import os


def _to_bool_local(value: Any, default: bool = False) -> bool:
    """Lightweight bool parser used only for strict-mode wiring.

    This intentionally mirrors the behaviour of the more general
    ``_to_bool`` helper in processor.py but is kept local to avoid
    import cycles.
    """

    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.lower() in ("true", "1", "yes", "on", "t", "y")
    if isinstance(value, (int, float)):
        return bool(value)
    return default


def configure_strict_science(merged_options: Dict[str, Any]) -> Tuple[Dict[str, Any], bool]:
    """Derive and apply the strict-science flag.

    Priority order (highest to lowest):
      1. ``merged_options['strict_science']`` if explicitly provided.
      2. ``SARDINE_STRICT`` environment variable (legacy integrations).
      3. Default: True (strict scientific mode).

    When strict science is enabled this helper:
      * Ensures ``SARDINE_STRICT=1`` is visible to the Rust core.
      * Forces validation and orbit handling into fail-fast mode by
        setting sensible defaults in the returned options dict.

    The returned ``strict_science`` boolean is intended to be stored on
    the BackscatterProcessor instance for downstream stages.
    """

    # Make a shallow copy so callers can safely reuse their original dict
    # without accidental mutation.
    opts = dict(merged_options) if not isinstance(merged_options, dict) else merged_options

    raw_opt_strict = opts.get("strict_science", None)
    if raw_opt_strict is not None:
        strict_science = _to_bool_local(raw_opt_strict, default=True)
    else:
        env_val = os.environ.get("SARDINE_STRICT")
        if env_val is not None:
            strict_science = _to_bool_local(env_val, default=True)
        else:
            strict_science = True

    if strict_science:
        # Ensure the Rust core and any helper scripts see strict mode.
        os.environ["SARDINE_STRICT"] = "1"
        opts["strict_science"] = True
        # Enforce fail-fast behaviour for validation and orbit.
        opts.setdefault("validation_enabled", True)
        opts.setdefault("validation_fail_fast", True)
        opts.setdefault("strict_orbit_mode", True)

    return opts, strict_science
