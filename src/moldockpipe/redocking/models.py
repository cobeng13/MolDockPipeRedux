from __future__ import annotations

from dataclasses import asdict, dataclass, field

from enum import Enum

from ..parallel import logical_cpu_count


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
    cpu_count: int = field(default_factory=logical_cpu_count)

    def as_dict(self) -> dict[str, int | float]:
        return asdict(self)
