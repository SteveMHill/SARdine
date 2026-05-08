"""Common processor utilities shared across SAR pipelines."""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, Optional


class BaseProcessor:
    """Reusable helper mixin for SAR processing pipelines."""

    def __init__(self, output_dir: Path | str, options: Optional[Dict[str, Any]] = None) -> None:
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.options: Dict[str, Any] = options.copy() if options else {}
        self.start_time = time.time()
        self.processing_log: list[Dict[str, Any]] = []

    # NOTE: Kept lightweight so pipelines can override formatting if needed.
    def announce_step(self, step_num: int, step_name: str, details: str = "", status: str = "start") -> None:
        """Emit a console banner announcing the upcoming processing step."""

        status_normalized = (status or "start").lower()
        if status_normalized == "skip":
            icon, label = "ℹ️", "SKIP"
        elif status_normalized == "resume":
            icon, label = "🔄", "RESUME"
        else:
            icon, label = "➡️", "START"

        print(f"\n{icon}  {label} STEP {step_num}: {step_name}")
        if details:
            print(f"   {details}")

    def log_step(
        self,
        step_num: int,
        step_name: str,
        status: str,
        details: str = "",
        duration: Optional[float] = None,
    ) -> None:
        """Persist step timing/metadata for progress tracking."""

        log_entry: Dict[str, Any] = {
            "step": step_num,
            "name": step_name,
            "status": status,
            "details": details,
            "duration": duration,
            "timestamp": time.time(),
        }
        self.processing_log.append(log_entry)

        try:
            progress_path = self.output_dir / "processing_log_progress.json"
            with open(progress_path, "w", encoding="utf-8") as progress_file:
                json.dump(self.processing_log, progress_file, indent=2)

            stream_path = self.output_dir / "processing_log.jsonl"
            with open(stream_path, "a", encoding="utf-8") as stream_file:
                stream_file.write(json.dumps(log_entry) + "\n")
        except Exception as exc:  # pragma: no cover - logging failures should never abort processing
            print(f"   ⚠️  Failed to persist progress log: {exc}")

        duration_suffix = f" ({duration:.1f}s)" if duration else ""
        status_icon_map = {
            "success": "✅",
            "error": "❌",
            "warning": "⚠️",
            "skipped": "⏭️",
            "fallback": "ℹ️",
            "info": "ℹ️",
        }
        status_icon = status_icon_map.get(status, "🔄")
        print(f"   {status_icon} STEP {step_num}: {step_name}{duration_suffix}")
        if details:
            print(f"      {details}")

