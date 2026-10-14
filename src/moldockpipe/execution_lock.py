from __future__ import annotations

import os
from contextlib import contextmanager
from functools import wraps
from pathlib import Path
from threading import local


_owned = local()


@contextmanager
def project_execution_lock(root: Path):
    """One pipeline writer per project, across threads and processes.

    The OS releases the lock on process exit, including crashes. Nested stages
    on the owning thread may reuse it; another thread must acquire its own lock.
    """
    root = Path(root).resolve()
    key = os.path.normcase(str(root))
    owned = getattr(_owned, "projects", None)
    if owned is None:
        owned = _owned.projects = set()
    if key in owned:
        yield
        return
    handle = (root / ".pipeline.lock").open("a+b")
    try:
        try:
            if os.name == "nt":
                import msvcrt
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise RuntimeError("This project is already being processed in another window or terminal. "
                               "Wait for that run to finish or stop it first.") from exc
        owned.add(key)
        try:
            yield
        finally:
            owned.remove(key)
    finally:
        # Closing the descriptor releases its OS lock, including on exceptions.
        handle.close()


def exclusive_project(method):
    @wraps(method)
    def locked(self, *args, **kwargs):
        root = self.root if hasattr(self, "root") else self.repository.root
        with project_execution_lock(root):
            return method(self, *args, **kwargs)
    return locked
