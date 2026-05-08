"""Lightweight validation gates for SARdine processing stages.

Provides structured recording of per-stage checks with optional fail-fast
behaviour. These gates are intentionally simple and low-overhead so they can
run inside Python control-flow without affecting Rust hot loops.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, List, Optional
from numbers import Number


@dataclass
class ValidationRecord:
    stage: str
    metric: str
    value: Any
    expected: Any = None
    tolerance: Optional[float] = None
    passed: bool = True
    severity: str = "info"  # info|warning|error
    message: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


class ValidationError(RuntimeError):
    """Raised when a validation gate fails in fail-fast mode."""


class ValidationGates:
    """Collects validation checks and optionally enforces fail-fast semantics."""

    def __init__(self, enabled: bool = True, fail_fast: bool = True) -> None:
        self.enabled = bool(enabled)
        self.fail_fast = bool(fail_fast)
        self.records: List[ValidationRecord] = []

    def record(
        self,
        stage: str,
        metric: str,
        *,
        value: Any,
        expected: Any = None,
        tolerance: Optional[float] = None,
        passed: bool = True,
        severity: str = "info",
        message: str = "",
    ) -> ValidationRecord:
        record = ValidationRecord(
            stage=stage,
            metric=metric,
            value=value,
            expected=expected,
            tolerance=tolerance,
            passed=bool(passed),
            severity=severity,
            message=message,
        )
        if self.enabled:
            self.records.append(record)
        return record

    def _should_raise(self, passed: bool, severity: str) -> bool:
        if not self.enabled:
            return False
        if passed:
            return False
        if severity != "error":
            return False
        return self.fail_fast

    def require(
        self,
        condition: bool,
        stage: str,
        metric: str,
        *,
        value: Any = None,
        expected: Any = None,
        tolerance: Optional[float] = None,
        severity: str = "error",
        message: str = "",
    ) -> bool:
        passed = bool(condition)
        record = self.record(
            stage,
            metric,
            value=value,
            expected=expected,
            tolerance=tolerance,
            passed=passed,
            severity=severity,
            message=message,
        )
        if self._should_raise(passed, severity):
            raise ValidationError(message or f"Validation failed: {stage}.{metric} -> {value}")
        return passed

    def expect_within(
        self,
        stage: str,
        metric: str,
        *,
        value: Any,
        expected: Any,
        tolerance: float,
        severity: str = "error",
        message: str = "",
    ) -> bool:
        passed = False
        if isinstance(value, Number) and isinstance(expected, Number):
            passed = abs(float(value) - float(expected)) <= float(tolerance)
        record_msg = message or ""
        self.record(
            stage,
            metric,
            value=value,
            expected=expected,
            tolerance=tolerance,
            passed=passed,
            severity=severity,
            message=record_msg,
        )
        if self._should_raise(passed, severity):
            raise ValidationError(record_msg or f"{stage}.{metric} outside tolerance: {value} vs {expected} ±{tolerance}")
        return passed

    def summary(self) -> dict:
        total = len(self.records)
        failed = sum(1 for r in self.records if not r.passed)
        warnings = sum(1 for r in self.records if r.severity == "warning")
        return {
            "total_checks": total,
            "failed": failed,
            "warnings": warnings,
            "passed": total - failed,
            "fail_fast": self.fail_fast,
            "enabled": self.enabled,
        }

    def save_json(self, path: str | Path) -> Path:
        output_path = Path(path)
        if not self.enabled:
            return output_path
        payload = {
            "summary": self.summary(),
            "records": [r.to_dict() for r in self.records],
        }
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        return output_path