from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path


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
