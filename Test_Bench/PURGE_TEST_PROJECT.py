#!/usr/bin/env python3
"""Fast, confirmation-protected cleanup for the disposable Test_Bench folder.

Usage: python Test_Bench/PURGE_TEST_PROJECT.py
       python Test_Bench/PURGE_TEST_PROJECT.py <project-folder>
The source CSV/receptor/config files are retained. Generated artifacts, logs,
exports, and SQLite workflow rows are removed so a one-compound test can be
repeated quickly. No third-party packages are required.
"""

from __future__ import annotations

import argparse
import shutil
import sqlite3
from pathlib import Path


GENERATED_DIRS = ("artifacts", "logs", "exports", "output", "3D_Structures", "prepared_ligands", "results", "state")


def purge(root: Path, *, confirmed: bool = False) -> None:
    root = root.resolve()
    if not confirmed:
        answer = input(f"Purge generated test data under {root}? Type PURGE to continue: ").strip()
        if answer != "PURGE":
            print("Cancelled.")
            return
    for name in GENERATED_DIRS:
        target = root / name
        if target.exists():
            shutil.rmtree(target)
        target.mkdir(parents=True, exist_ok=True)
    database = root / "project.sqlite"
    if database.exists():
        with sqlite3.connect(database) as conn:
            conn.execute("PRAGMA foreign_keys = OFF")
        for table in ("docking_poses", "docking_runs", "conformers", "molecular_states", "screening_results", "stage_runs", "artifacts", "manual_reviews", "parent_ligands"):
            conn.execute(f"DELETE FROM {table}")
        conn.commit()
    print(f"Purged generated test data; preserved {root / 'inputs'} and project settings.")


def main() -> int:
    parser = argparse.ArgumentParser()
    default_project = Path(__file__).resolve().parent
    parser.add_argument("project", type=Path, nargs="?", default=default_project,
                        help=f"Project to purge (default: {default_project})")
    parser.add_argument("--yes", action="store_true", help="Skip the confirmation prompt")
    args = parser.parse_args()
    purge(args.project, confirmed=args.yes)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
