"""Isolated dashboard scale probe. Never opens or changes a user project.

Run with the application's Python and QT_QPA_PLATFORM=offscreen.
Uses synthetic SQLite records; does not run chemistry or Vina.
"""
from __future__ import annotations

import argparse
import ctypes
import json
import os
from pathlib import Path
import statistics
import tempfile
import time
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt6.QtCore import QTimer
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication
from moldockpipe.project import DEFAULT_VINA_PROFILE, ProjectRepository
from moldockpipe.ui.main_window import MainWindow


def memory_mib():
    if os.name != "nt":
        return None
    from ctypes import wintypes

    class Counters(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD)] + [
            (name, ctypes.c_size_t) for name in (
                "PeakWorkingSetSize", "WorkingSetSize", "QuotaPeakPagedPoolUsage",
                "QuotaPagedPoolUsage", "QuotaPeakNonPagedPoolUsage", "QuotaNonPagedPoolUsage",
                "PagefileUsage", "PeakPagefileUsage", "PrivateUsage")]

    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.GetCurrentProcess.restype = wintypes.HANDLE
    counters = Counters()
    counters.cb = ctypes.sizeof(counters)
    psapi = ctypes.WinDLL("psapi", use_last_error=True)
    psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.c_void_p, wintypes.DWORD]
    if not psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(), ctypes.byref(counters), counters.cb):
        raise ctypes.WinError(ctypes.get_last_error())
    return {"working_set": round(counters.WorkingSetSize / 2**20, 2),
            "private": round(counters.PrivateUsage / 2**20, 2)}


def seed(repo, parents, states_per_parent):
    repo.save_receptor_profiles([{**DEFAULT_VINA_PROFILE, "id": "rec"}])
    with repo.connection() as conn:
        conn.executemany("INSERT INTO parent_ligands(parent_id,source_smiles,parse_status,created_at) VALUES (?,'CCO','completed','now')",
                         ((f"p{i}",) for i in range(parents)))
        conn.executemany("INSERT INTO screening_results(parent_id,fingerprint,status,descriptors_json,rules_json,decision,created_at) VALUES (?,'fp','completed','{}','{}','pass','now')",
                         ((f"p{i}",) for i in range(parents)))
        conn.executemany("INSERT INTO stage_runs VALUES (?,'molscrub','parent',?,'fp','completed','now','now',NULL)",
                         ((f"history{i}", f"p{i}") for i in range(parents)))
        conn.executemany("INSERT INTO molecular_states(state_id,parent_id,state_smiles,state_isomeric_smiles,formal_charge,state_structure_hash,generation_fingerprint,status,created_at) VALUES (?,?,'CCO','CCO',0,?,'fp','pdbqt_ready','now')",
                         ((f"s{i}-{j}", f"p{i}", f"h{j}") for i in range(parents) for j in range(states_per_parent)))
        conn.executemany("INSERT INTO conformers VALUES (?,?,1,NULL,'pdbqt_ready','now')",
                         ((f"c{i}-{j}", f"s{i}-{j}") for i in range(parents) for j in range(states_per_parent)))
        conn.executemany("INSERT INTO docking_runs(run_id,state_id,receptor_hash,settings_fingerprint,status,is_current,receptor_profile_id) VALUES (?,?,'hash','fp','completed',1,'rec')",
                         ((f"r{i}-{j}", f"s{i}-{j}") for i in range(parents) for j in range(states_per_parent)))
        conn.executemany("INSERT INTO docking_poses(pose_id,run_id,mode_index,affinity) VALUES (?,?,?,-7.0)",
                         ((f"pose{i}-{j}-{k}", f"r{i}-{j}", k) for i in range(parents) for j in range(states_per_parent) for k in range(1, 10)))


def run_case(app, parents, states_per_parent, repeats):
    with tempfile.TemporaryDirectory(prefix="moldock-dashboard-") as temporary:
        repo = ProjectRepository.create(Path(temporary) / "project")
        seed(repo, parents, states_per_parent)
        window = MainWindow()
        window.repo = repo
        window.stage_running = True
        window.show()
        app.processEvents()
        times, gaps, request_times = [], [], []
        last = time.perf_counter()

        def tick():
            nonlocal last
            now = time.perf_counter()
            gaps.append(now - last)
            last = now

        timer = QTimer()
        timer.timeout.connect(tick)
        timer.start(10)
        before = memory_mib()
        for _ in range(repeats):
            started = time.perf_counter()
            window._refresh_live_dashboard()
            request_times.append(time.perf_counter() - started)
            deadline = started + 120
            while window._dashboard_worker is not None:
                app.processEvents()
                QTest.qWait(2)
                if time.perf_counter() > deadline:
                    raise TimeoutError("Dashboard benchmark exceeded 120s")
            times.append(window._snapshot_seconds)
            # Force the next measured sample without including the intentional
            # completion-based one-second cooldown in collection timings.
            window._dashboard_timer.stop()
        after = memory_mib()
        timer.stop()
        # Bound the legacy query rather than hanging the investigation.
        with repo.connection() as conn:
            conn.execute("DROP INDEX idx_stage_runs_entity")
        deadline = time.perf_counter() + 3
        started = time.perf_counter()
        with repo.connection() as conn:
            conn.set_progress_handler(lambda: int(time.perf_counter() > deadline), 10000)
            try:
                conn.execute("""SELECT p.parent_id FROM parent_ligands p LEFT JOIN stage_runs r ON r.rowid=(
                    SELECT MAX(latest.rowid) FROM stage_runs latest WHERE latest.stage_name='molscrub'
                    AND latest.entity_type='parent' AND latest.entity_id=p.parent_id)
                    WHERE p.active=1 AND r.status='failed'""").fetchall()
                legacy = time.perf_counter() - started
            except Exception as exc:
                if "interrupted" not in str(exc):
                    raise
                legacy = ">3s (interrupted)"
            finally:
                conn.set_progress_handler(None, 0)
        window.stage_running = False
        window.close()
        app.processEvents()
        return {"parents": parents, "states_per_parent": states_per_parent,
                "poses": parents * states_per_parent * 9, "repeats": repeats,
                "refresh_seconds_median": round(statistics.median(times), 4),
                "refresh_seconds_max": round(max(times), 4),
                "gui_request_seconds_max": round(max(request_times), 4),
                "qt_heartbeat_max_gap_seconds": round(max(gaps, default=0), 4),
                "memory_before_mib": before, "memory_after_mib": after,
                "legacy_unindexed_failure_query_seconds": legacy}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=10)
    parser.add_argument("--case", type=int, nargs=2, action="append", metavar=("PARENTS", "STATES"))
    args = parser.parse_args()
    app = QApplication([])
    results = []
    for parents, states in (args.case or ((1000, 1), (5000, 1), (20000, 1), (20000, 8))):
        result = run_case(app, parents, states, args.repeats)
        results.append(result)
        print(json.dumps(result), flush=True)
        args.output.write_text(json.dumps({"python": sys.version, "cpu_count": os.cpu_count(),
                                          "cases": results}, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
