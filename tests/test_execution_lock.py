import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from moldockpipe.execution_lock import project_execution_lock


def test_project_lock_is_reentrant_but_excludes_other_threads(tmp_path):
    with project_execution_lock(tmp_path):
        with project_execution_lock(tmp_path):
            with ThreadPoolExecutor(max_workers=1) as pool:
                future = pool.submit(lambda: _acquire(tmp_path))
                with pytest.raises(RuntimeError, match="already being processed"):
                    future.result(timeout=5)
    _acquire(tmp_path)


def _acquire(root):
    with project_execution_lock(root):
        return True


def _child(root, code):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1] / "src")
    return subprocess.run([sys.executable, "-B", "-c", code, str(root)],
                          env=env, capture_output=True, text=True, timeout=15, check=True)


def test_project_lock_excludes_other_processes_and_releases_after_crash(tmp_path):
    attempt = """import sys
from pathlib import Path
from moldockpipe.execution_lock import project_execution_lock
try:
    with project_execution_lock(Path(sys.argv[1])):
        print('acquired')
except RuntimeError:
    print('busy')
"""
    with project_execution_lock(tmp_path):
        assert _child(tmp_path, attempt).stdout.strip() == "busy"
    assert _child(tmp_path, attempt).stdout.strip() == "acquired"
    crash = """import os,sys
from pathlib import Path
from moldockpipe.execution_lock import project_execution_lock
with project_execution_lock(Path(sys.argv[1])):
    os._exit(0)
"""
    _child(tmp_path, crash)
    _acquire(tmp_path)


def test_different_projects_can_run_concurrently(tmp_path):
    first, second = tmp_path / "first", tmp_path / "second"
    first.mkdir()
    second.mkdir()
    with project_execution_lock(first):
        with ThreadPoolExecutor(max_workers=1) as pool:
            assert pool.submit(_acquire, second).result(timeout=5)


def test_project_lock_releases_after_failure(tmp_path):
    with pytest.raises(ValueError):
        with project_execution_lock(tmp_path):
            raise ValueError("failed preparation")
    _acquire(tmp_path)
