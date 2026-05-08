"""Lightweight orchestration primitives for SARdine processing pipelines."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, List, MutableMapping, Optional, Sequence
import logging
import time

_log = logging.getLogger(__name__)


@dataclass
class PipelineContext:
    """Shared execution context passed between processing stages."""

    processor: Any
    state: MutableMapping[str, Any] = field(default_factory=dict)
    artifacts: MutableMapping[str, Any] = field(default_factory=dict)

    def set_artifact(self, key: str, value: Any) -> None:
        self.artifacts[key] = value

    def get_artifact(self, key: str, default: Any = None) -> Any:
        return self.artifacts.get(key, default)


class PipelineStage:
    """Base stage abstraction. Subclasses override ``run``."""

    name: str = "stage"
    description: str = ""

    def run(self, context: PipelineContext) -> None:  # pragma: no cover - interface only
        raise NotImplementedError


class FunctionStage(PipelineStage):
    """Stage adapter around a simple callable."""

    def __init__(self, name: str, func: Callable[[PipelineContext], None], description: str = "") -> None:
        self.name = name
        self.description = description
        self._func = func

    def run(self, context: PipelineContext) -> None:
        self._func(context)


class ProcessingPipeline:
    """Sequential pipeline executor for reusable stage composition."""

    def __init__(self, stages: Sequence[PipelineStage]) -> None:
        self._stages: List[PipelineStage] = list(stages)

    @property
    def stages(self) -> Sequence[PipelineStage]:
        return tuple(self._stages)

    def run(self, context: PipelineContext) -> None:
        processor = getattr(context, "processor", None)
        for index, stage in enumerate(self._stages, start=1):
            context.state["stage_index"] = index
            context.state["stage_name"] = stage.name
            context.state["stage_description"] = stage.description
            start = time.perf_counter()
            stage.run(context)
            duration = time.perf_counter() - start

            if processor is not None:
                # Record timing for downstream summaries
                record = getattr(processor, "_record_stage_timing", None)
                if callable(record):
                    try:
                        record(stage.name, duration)
                    except Exception:
                        _log.warning("Failed to record stage timing for '%s'", stage.name, exc_info=True)
                # Optional verbose trace for live runs
                logger = getattr(processor, "_log", None)
                if callable(logger):
                    try:
                        logger(f"\u23f1\ufe0f  Stage {index} ({stage.name}) finished in {duration:.2f}s", level=2)
                    except Exception:
                        _log.warning("Failed to log stage completion for '%s'", stage.name, exc_info=True)
                # Run post-stage validation if available
                validate = getattr(processor, "_validate_stage_output", None)
                if callable(validate):
                    try:
                        validate(stage.name, context)
                    except Exception as e:
                        # Don't fail pipeline on validation errors, just log
                        if logger and callable(logger):
                            logger(f"⚠️  Validation warning after {stage.name}: {e}", level=1)
                        pass


__all__ = [
    "PipelineContext",
    "PipelineStage",
    "FunctionStage",
    "ProcessingPipeline",
]
