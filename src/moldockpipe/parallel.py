from __future__ import annotations

import os
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from typing import Callable, Iterable, TypeVar


T = TypeVar("T")
R = TypeVar("R")


def logical_cpu_count() -> int:
    return max(1, os.cpu_count() or 1)


def worker_count(settings: dict[str, object], jobs: int) -> int:
    """Zero means automatic; explicit positive limits remain portable."""
    requested = int(settings.get("workers", 0))
    if requested < 0:
        raise ValueError("Parallel workers must be zero (automatic) or a positive number")
    return max(1, min(requested or logical_cpu_count(), max(1, jobs)))


def bounded_futures(function: Callable[[T], R], items: Iterable[T], workers: int,
                    *, cancelled: Callable[[], bool] | None = None, name: str = "preparation",
                    on_discard: Callable[[T], None] | None = None):
    """Stop submitting on cancellation and drain active work safely.

    The caller consumes results so database writes and progress notifications
    stay on the coordinating pipeline thread. Only one wave is queued.
    """
    iterator = iter(items)
    exhausted = False
    with ThreadPoolExecutor(max_workers=workers, thread_name_prefix=name) as pool:
        pending = {}
        while pending or not exhausted:
            while len(pending) < workers and not exhausted and not (cancelled and cancelled()):
                try:
                    item = next(iterator)
                except StopIteration:
                    exhausted = True
                    break
                # Selecting a job can do I/O and emit UI events. Stop may arrive
                # during that work, after the loop's first cancellation check.
                if cancelled and cancelled():
                    if on_discard:
                        on_discard(item)
                    exhausted = True
                    break
                pending[pool.submit(function, item)] = item
            if not pending:
                break
            done, _ = wait(pending, timeout=0.1, return_when=FIRST_COMPLETED)
            for future in done:
                yield pending.pop(future), future
    if cancelled and cancelled():
        raise InterruptedError("Pipeline stopped safely")
