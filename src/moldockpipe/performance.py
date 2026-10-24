"""On-demand, local performance diagnostic snapshot without chemical structures."""
import json
import hashlib
import os
import platform
import sys
import time
from pathlib import Path
from datetime import datetime, timezone
from .resources import available_memory_bytes, process_tree_snapshot
from . import __version__


def save_performance_report(repo, path, gui_metrics):
    started = time.monotonic()
    with repo.connection() as conn:
        counts = {table: conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
                  for table in ('parent_ligands', 'molecular_states', 'conformers',
                                'docking_runs', 'docking_poses', 'stage_runs')}
        indexes = {table: [row[1] for row in conn.execute(f"PRAGMA index_list('{table}')")]
                   for table in ('stage_runs', 'molecular_states', 'conformers', 'docking_runs', 'docking_poses')}
    source = Path(__file__).parent
    report = {"captured_at": datetime.now(timezone.utc).isoformat(),
              "application_version": __version__,
              "source_sha256": {name: hashlib.sha256((source / name).read_bytes()).hexdigest()
                                for name in ('dashboard.py', 'parallel.py', 'ui/main_window.py', 'ui/workers.py')},
              "python": sys.executable, "python_version": sys.version,
              "platform": platform.platform(), "logical_cpus": os.cpu_count(),
              "project": str(repo.root), "source": str(__file__),
              "available_memory_bytes": available_memory_bytes(),
              "processes": process_tree_snapshot(), "gui": dict(gui_metrics),
              "record_counts": counts, "indexes": indexes,
              "collection_seconds": time.monotonic() - started,
              "note": "One local snapshot; CPU times are cumulative, not utilization percentages. No ligand structures included."}
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return path
