from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum


class RedockingStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    INTERRUPTED = "interrupted"
    FAILED = "failed"
    CANCELLED = "cancelled"
    ARTIFACTS_READY = "ARTIFACTS_READY"


@dataclass(frozen=True)
class RedockingSettings:
    exhaustiveness: int = 32
    num_modes: int = 20
    energy_range: float = 5.0
    seed: int = 123456
    cpu_count: int = 1

    def as_dict(self) -> dict[str, int | float]:
        return asdict(self)
