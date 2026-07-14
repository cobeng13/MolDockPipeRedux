from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any, Callable


class StageStatus(StrEnum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"
    INTERRUPTED = "interrupted"
    SKIPPED = "skipped"


class ScreeningPolicy(StrEnum):
    ANNOTATE_ONLY = "annotate_only"
    EXCLUDE_FAILING_ALL = "exclude_failing_all"
    EXCLUDE_FAILING_ANY = "exclude_failing_any"
    MANUAL_REVIEW = "manual_review"


@dataclass(frozen=True)
class Artifact:
    path: Path
    sha256: str
    size: int
    kind: str


@dataclass(frozen=True)
class ProgressEvent:
    event: str
    stage: str
    item_id: str | None = None
    index: int = 0
    total: int = 0
    succeeded: int = 0
    skipped: int = 0
    failed: int = 0
    message: str = ""
    receptor_profile_id: str | None = None
    receptor_profile_name: str | None = None


ProgressCallback = Callable[[ProgressEvent], None]
