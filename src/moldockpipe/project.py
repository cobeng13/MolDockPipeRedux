from __future__ import annotations

import csv
import json
import platform
import sqlite3
import sys
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

import yaml

from .fingerprints import file_sha256
from .models import ScreeningPolicy, StageStatus


SCHEMA_VERSION = 1


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


class ProjectRepository:
    """Single-writer access point for a portable docking project."""

    def __init__(self, root: Path):
        self.root = root.resolve()
        self.database_path = self.root / "project.sqlite"

    @classmethod
    def create(cls, root: Path, name: str | None = None) -> "ProjectRepository":
        root = root.resolve()
        root.mkdir(parents=True, exist_ok=True)
        for relative in ("inputs", "artifacts/sdf", "artifacts/pdbqt", "artifacts/docking", "logs", "exports", "tools/vina"):
            (root / relative).mkdir(parents=True, exist_ok=True)
        config = {
            "name": name or root.name,
            "schema_version": SCHEMA_VERSION,
            "screening": {"policy": ScreeningPolicy.ANNOTATE_ONLY.value, "lipinski": True, "veber": True, "egan": False, "ghose": False},
            "molscrub": {"ph": 7.4, "enumerate_states": True, "max_states": 32, "stereochemistry_policy": "warn_and_continue", "fragment_policy": "manual_review"},
            "vina": {"receptor": "inputs/receptor_prepared.pdbqt", "center_x": 0.0, "center_y": 0.0, "center_z": 0.0, "size_x": 20.0, "size_y": 20.0, "size_z": 20.0, "exhaustiveness": 8, "num_modes": 9, "energy_range": 3, "seed": 42, "cpu_count": 1},
        }
        (root / "project.yml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
        repo = cls(root)
        repo.migrate()
        repo.record_environment()
        return repo

    @contextmanager
    def connection(self) -> Iterator[sqlite3.Connection]:
        conn = sqlite3.connect(self.database_path)
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA foreign_keys = ON")
        conn.execute("PRAGMA journal_mode = WAL")
        try:
            yield conn
            conn.commit()
        except Exception:
            conn.rollback()
            raise
        finally:
            conn.close()

    def migrate(self) -> None:
        with self.connection() as conn:
            conn.executescript("""
            CREATE TABLE IF NOT EXISTS schema_migrations (version INTEGER PRIMARY KEY, applied_at TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS projects (project_id TEXT PRIMARY KEY, name TEXT NOT NULL, root_path TEXT NOT NULL, created_at TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS parent_ligands (
                parent_id TEXT PRIMARY KEY, source_id TEXT, source_smiles TEXT NOT NULL,
                canonical_source_smiles TEXT, parent_inchikey TEXT, notes TEXT, params_json TEXT,
                parse_status TEXT NOT NULL, parse_reason TEXT, created_at TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS screening_results (
                parent_id TEXT PRIMARY KEY REFERENCES parent_ligands(parent_id), fingerprint TEXT NOT NULL,
                status TEXT NOT NULL, descriptors_json TEXT NOT NULL, rules_json TEXT NOT NULL,
                decision TEXT NOT NULL, reason TEXT, created_at TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS molecular_states (
                state_id TEXT PRIMARY KEY, parent_id TEXT NOT NULL REFERENCES parent_ligands(parent_id),
                state_smiles TEXT NOT NULL, state_isomeric_smiles TEXT NOT NULL, state_inchikey TEXT,
                formal_charge INTEGER NOT NULL, tautomer_index INTEGER, protomer_index INTEGER,
                state_structure_hash TEXT NOT NULL, generation_fingerprint TEXT NOT NULL,
                status TEXT NOT NULL, reason TEXT, enumeration_truncated INTEGER NOT NULL DEFAULT 0,
                created_at TEXT NOT NULL, UNIQUE(parent_id, state_structure_hash, generation_fingerprint)
            );
            CREATE TABLE IF NOT EXISTS conformers (
                conformer_id TEXT PRIMARY KEY, state_id TEXT NOT NULL REFERENCES molecular_states(state_id),
                conformer_index INTEGER NOT NULL, sdf_artifact_id TEXT, status TEXT NOT NULL, created_at TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS stage_runs (
                stage_run_id TEXT PRIMARY KEY, stage_name TEXT NOT NULL, entity_type TEXT NOT NULL, entity_id TEXT NOT NULL,
                fingerprint TEXT NOT NULL, status TEXT NOT NULL, started_at TEXT, ended_at TEXT, reason TEXT
            );
            CREATE TABLE IF NOT EXISTS artifacts (
                artifact_id TEXT PRIMARY KEY, relative_path TEXT NOT NULL UNIQUE, sha256 TEXT NOT NULL, size_bytes INTEGER NOT NULL,
                artifact_type TEXT NOT NULL, created_stage TEXT NOT NULL, validation_status TEXT NOT NULL,
                active INTEGER NOT NULL DEFAULT 1, created_at TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS docking_runs (
                run_id TEXT PRIMARY KEY, state_id TEXT NOT NULL REFERENCES molecular_states(state_id), receptor_hash TEXT NOT NULL,
                settings_fingerprint TEXT NOT NULL, status TEXT NOT NULL, log_artifact_id TEXT REFERENCES artifacts(artifact_id),
                raw_output_artifact_id TEXT REFERENCES artifacts(artifact_id), command_json TEXT, return_code INTEGER,
                started_at TEXT, ended_at TEXT, is_current INTEGER NOT NULL DEFAULT 0
            );
            CREATE TABLE IF NOT EXISTS docking_poses (
                pose_id TEXT PRIMARY KEY, run_id TEXT NOT NULL REFERENCES docking_runs(run_id), mode_index INTEGER NOT NULL,
                affinity REAL NOT NULL, rmsd_lb REAL, rmsd_ub REAL, exported_sdf_record INTEGER,
                is_best_for_state INTEGER NOT NULL DEFAULT 0, UNIQUE(run_id, mode_index)
            );
            CREATE TABLE IF NOT EXISTS tool_installations (tool_name TEXT PRIMARY KEY, version TEXT, location TEXT, sha256 TEXT, recorded_at TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS project_settings (key TEXT PRIMARY KEY, value_json TEXT NOT NULL, updated_at TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS manual_reviews (review_id TEXT PRIMARY KEY, parent_id TEXT REFERENCES parent_ligands(parent_id), decision TEXT NOT NULL, note TEXT, created_at TEXT NOT NULL);
            """)
            conn.execute("INSERT OR IGNORE INTO schema_migrations(version, applied_at) VALUES (?, ?)", (SCHEMA_VERSION, utc_now()))
            conn.execute("INSERT OR IGNORE INTO projects(project_id, name, root_path, created_at) VALUES (?, ?, ?, ?)", (str(uuid.uuid4()), self.root.name, str(self.root), utc_now()))

    def record_environment(self) -> None:
        self.record_tool("python", platform.python_version(), sys.executable)
        self.set_setting("environment", {"os": platform.platform(), "architecture": platform.machine(), "python": sys.version})

    def record_tool(self, name: str, version: str | None, location: str | Path | None) -> None:
        path = Path(location) if location else None
        digest = file_sha256(path) if path and path.is_file() else None
        with self.connection() as conn:
            conn.execute("INSERT OR REPLACE INTO tool_installations VALUES (?, ?, ?, ?, ?)", (name, version, str(location or ""), digest, utc_now()))

    def set_setting(self, key: str, value: Any) -> None:
        with self.connection() as conn:
            conn.execute("INSERT OR REPLACE INTO project_settings VALUES (?, ?, ?)", (key, json.dumps(value, sort_keys=True), utc_now()))

    def get_settings(self) -> dict[str, Any]:
        path = self.root / "project.yml"
        return yaml.safe_load(path.read_text(encoding="utf-8")) if path.exists() else {}

    def import_ligands_csv(self, path: Path) -> int:
        count = 0
        with path.open(newline="", encoding="utf-8-sig") as handle, self.connection() as conn:
            for index, row in enumerate(csv.DictReader(handle), 1):
                smiles = (row.get("smiles") or "").strip()
                if not smiles:
                    continue
                parent_id = row.get("id", "").strip() or f"ligand_{index:05d}"
                conn.execute("""INSERT OR REPLACE INTO parent_ligands
                    (parent_id, source_id, source_smiles, notes, params_json, parse_status, created_at)
                    VALUES (?, ?, ?, ?, ?, ?, ?)""", (parent_id, row.get("id"), smiles, row.get("notes"), row.get("params_json"), StageStatus.PENDING, utc_now()))
                count += 1
        return count

    def clear_workflow_data(self) -> None:
        """Clear current ligand/stage data while retaining project settings."""
        with self.connection() as conn:
            conn.execute("PRAGMA foreign_keys = OFF")
            for table in ("docking_poses", "docking_runs", "conformers", "molecular_states", "screening_results", "stage_runs", "artifacts", "manual_reviews", "parent_ligands"):
                conn.execute(f"DELETE FROM {table}")
            conn.execute("PRAGMA foreign_keys = ON")

    def add_artifact(self, path: Path, artifact_type: str, stage: str, validation_status: str = "valid") -> str:
        path = path.resolve()
        artifact_id = str(uuid.uuid4())
        relative = path.relative_to(self.root).as_posix()
        with self.connection() as conn:
            conn.execute("INSERT INTO artifacts VALUES (?, ?, ?, ?, ?, ?, ?, 1, ?)", (artifact_id, relative, file_sha256(path), path.stat().st_size, artifact_type, stage, validation_status, utc_now()))
        return artifact_id

    def recover_interrupted_runs(self) -> int:
        with self.connection() as conn:
            result = conn.execute("UPDATE stage_runs SET status=?, ended_at=?, reason=? WHERE status=?", (StageStatus.INTERRUPTED, utc_now(), "Application restarted before completion", StageStatus.RUNNING)).rowcount
            conn.execute("UPDATE docking_runs SET status=?, ended_at=? WHERE status=?", (StageStatus.INTERRUPTED, utc_now(), StageStatus.RUNNING))
            return result
