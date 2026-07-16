from __future__ import annotations

import csv
import json
import platform
import shutil
import sqlite3
import sys
import uuid
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

import yaml

from .fingerprints import file_sha256, fingerprint
from .models import ScreeningPolicy, StageStatus


SCHEMA_VERSION = 4


DEFAULT_VINA_PROFILE: dict[str, Any] = {
    "id": "default",
    "name": "Default",
    "enabled": True,
    "archived": False,
    "receptor": "inputs/receptors/default/receptor.pdbqt",
    "center_x": 0.0,
    "center_y": 0.0,
    "center_z": 0.0,
    "size_x": 20.0,
    "size_y": 20.0,
    "size_z": 20.0,
    "exhaustiveness": 8,
    "num_modes": 9,
    "energy_range": 3,
    "seed": 42,
    "cpu_count": 1,
}


@dataclass(frozen=True)
class LigandSyncResult:
    added: int = 0
    unchanged: int = 0
    changed: int = 0
    archived: int = 0
    rejected: int = 0

    @property
    def reprocess_required(self) -> int:
        return self.added + self.changed


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
        for relative in ("inputs", "inputs/receptors/default", "artifacts/sdf", "artifacts/pdbqt", "artifacts/docking", "logs", "exports", "For_PostDocking", "tools/vina"):
            (root / relative).mkdir(parents=True, exist_ok=True)
        config = {
            "name": name or root.name,
            "schema_version": SCHEMA_VERSION,
            "screening": {"policy": ScreeningPolicy.ANNOTATE_ONLY.value, "lipinski": True, "veber": True, "egan": False, "ghose": False, "boiled_egg": False},
            "molscrub": {"ph": 7.4, "enumerate_states": True, "max_states": 32, "stereochemistry_policy": "warn_and_continue", "fragment_policy": "manual_review"},
            "meeko": {"workers": 4},
            "postdock": {"mode": "split_and_sdf", "poses_per_compound": 3, "selected_parents": []},
            "vina": {"profiles": [dict(DEFAULT_VINA_PROFILE)]},
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
                parse_status TEXT NOT NULL, parse_reason TEXT, created_at TEXT NOT NULL,
                active INTEGER NOT NULL DEFAULT 1, input_fingerprint TEXT, archived_at TEXT
            );
            CREATE TABLE IF NOT EXISTS screening_results (
                parent_id TEXT PRIMARY KEY REFERENCES parent_ligands(parent_id), fingerprint TEXT NOT NULL,
                status TEXT NOT NULL, descriptors_json TEXT NOT NULL, rules_json TEXT NOT NULL,
                decision TEXT NOT NULL, reason TEXT, created_at TEXT NOT NULL,
                active INTEGER NOT NULL DEFAULT 1
            );
            CREATE TABLE IF NOT EXISTS molecular_states (
                state_id TEXT PRIMARY KEY, parent_id TEXT NOT NULL REFERENCES parent_ligands(parent_id),
                state_smiles TEXT NOT NULL, state_isomeric_smiles TEXT NOT NULL, state_inchikey TEXT,
                formal_charge INTEGER NOT NULL, tautomer_index INTEGER, protomer_index INTEGER,
                state_structure_hash TEXT NOT NULL, generation_fingerprint TEXT NOT NULL,
                status TEXT NOT NULL, reason TEXT, enumeration_truncated INTEGER NOT NULL DEFAULT 0,
                created_at TEXT NOT NULL, active INTEGER NOT NULL DEFAULT 1,
                UNIQUE(parent_id, state_structure_hash, generation_fingerprint)
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
            CREATE TABLE IF NOT EXISTS redocking_runs (
                run_id TEXT PRIMARY KEY, receptor_profile_id TEXT NOT NULL, status TEXT NOT NULL,
                reference_ligand_id TEXT NOT NULL, receptor_path TEXT NOT NULL, receptor_sha256 TEXT NOT NULL,
                reference_sdf_path TEXT NOT NULL, reference_mol2_path TEXT NOT NULL, reference_sha256 TEXT NOT NULL,
                prepared_ligand_path TEXT, prepared_ligand_sha256 TEXT, settings_json TEXT NOT NULL,
                fingerprints_json TEXT NOT NULL DEFAULT '{}', meeko_version TEXT, vina_version TEXT,
                current_stage TEXT, started_at TEXT, finished_at TEXT, interrupted_at TEXT,
                error_stage TEXT, error_message TEXT, cancel_requested INTEGER NOT NULL DEFAULT 0
            );
            CREATE TABLE IF NOT EXISTS redocking_stages (
                run_id TEXT NOT NULL REFERENCES redocking_runs(run_id), stage_name TEXT NOT NULL,
                fingerprint TEXT NOT NULL, status TEXT NOT NULL, artifacts_json TEXT NOT NULL DEFAULT '{}',
                started_at TEXT, finished_at TEXT, error_message TEXT, PRIMARY KEY(run_id, stage_name)
            );
            CREATE TABLE IF NOT EXISTS redocking_poses (
                pose_id TEXT PRIMARY KEY, run_id TEXT NOT NULL REFERENCES redocking_runs(run_id), pose_rank INTEGER NOT NULL,
                affinity REAL, sdf_path TEXT NOT NULL, mol2_path TEXT NOT NULL, pdbqt_source_path TEXT NOT NULL,
                artifact_sha256 TEXT NOT NULL, metadata_json TEXT NOT NULL DEFAULT '{}', UNIQUE(run_id, pose_rank)
            );
            CREATE INDEX IF NOT EXISTS idx_redocking_profile_status ON redocking_runs(receptor_profile_id, status);
            """)
            self._ensure_column(conn, "parent_ligands", "active", "INTEGER NOT NULL DEFAULT 1")
            self._ensure_column(conn, "parent_ligands", "input_fingerprint", "TEXT")
            self._ensure_column(conn, "parent_ligands", "archived_at", "TEXT")
            self._ensure_column(conn, "screening_results", "active", "INTEGER NOT NULL DEFAULT 1")
            self._ensure_column(conn, "molecular_states", "active", "INTEGER NOT NULL DEFAULT 1")
            self._ensure_column(conn, "docking_runs", "receptor_profile_id", "TEXT NOT NULL DEFAULT 'default'")
            self._ensure_column(conn, "docking_runs", "receptor_profile_name", "TEXT NOT NULL DEFAULT 'Default'")
            self._ensure_column(conn, "docking_runs", "reason", "TEXT")
            conn.execute("""WITH ranked AS (
                    SELECT run_id, ROW_NUMBER() OVER (
                        PARTITION BY state_id, receptor_profile_id
                        ORDER BY CASE WHEN status='completed' THEN 0 ELSE 1 END,
                                 COALESCE(ended_at, started_at, '') DESC, run_id DESC) rank
                    FROM docking_runs WHERE is_current=1)
                UPDATE docking_runs SET is_current=0 WHERE run_id IN (SELECT run_id FROM ranked WHERE rank>1)""")
            conn.execute("CREATE INDEX IF NOT EXISTS idx_docking_runs_profile_state ON docking_runs(receptor_profile_id, state_id, is_current)")
            for row in conn.execute("SELECT parent_id, source_smiles, notes, params_json FROM parent_ligands WHERE input_fingerprint IS NULL OR input_fingerprint='' ").fetchall():
                conn.execute("UPDATE parent_ligands SET input_fingerprint=? WHERE parent_id=?", (self._input_fingerprint(row["source_smiles"], row["notes"], row["params_json"]), row["parent_id"]))
            conn.execute("INSERT OR IGNORE INTO schema_migrations(version, applied_at) VALUES (?, ?)", (SCHEMA_VERSION, utc_now()))
            conn.execute("INSERT OR IGNORE INTO projects(project_id, name, root_path, created_at) VALUES (?, ?, ?, ?)", (str(uuid.uuid4()), self.root.name, str(self.root), utc_now()))
        self._migrate_vina_profiles()

    def _migrate_vina_profiles(self) -> None:
        """Convert the legacy flat Vina block into a portable receptor profile."""
        path = self.root / "project.yml"
        if not path.exists():
            return
        config = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        vina = dict(config.get("vina", {}) or {})
        if isinstance(vina.get("profiles"), list):
            if config.get("schema_version") != SCHEMA_VERSION:
                config["schema_version"] = SCHEMA_VERSION
                path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
            return
        profile = dict(DEFAULT_VINA_PROFILE)
        legacy_receptor = str(vina.get("receptor", "inputs/receptor_prepared.pdbqt"))
        for key in ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z", "exhaustiveness", "num_modes", "energy_range", "seed", "cpu_count"):
            if key in vina:
                profile[key] = vina[key]
        source = self.root / legacy_receptor
        destination = self.root / str(profile["receptor"])
        if source.is_file() and source.resolve() != destination.resolve():
            destination.parent.mkdir(parents=True, exist_ok=True)
            if not destination.exists():
                shutil.copy2(source, destination)
        config["schema_version"] = SCHEMA_VERSION
        config["vina"] = {"profiles": [profile]}
        path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    def get_receptor_profiles(self, *, include_archived: bool = False) -> list[dict[str, Any]]:
        profiles = self.get_settings().get("vina", {}).get("profiles", [])
        if not isinstance(profiles, list):
            return []
        return [dict(profile) for profile in profiles if isinstance(profile, dict) and (include_archived or not profile.get("archived", False))]

    def save_receptor_profiles(self, profiles: list[dict[str, Any]]) -> None:
        path = self.root / "project.yml"
        config = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        config["vina"] = {"profiles": profiles}
        path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")

    @staticmethod
    def _ensure_column(conn: sqlite3.Connection, table: str, column: str, definition: str) -> None:
        columns = {row[1] for row in conn.execute(f"PRAGMA table_info({table})").fetchall()}
        if column not in columns:
            conn.execute(f"ALTER TABLE {table} ADD COLUMN {column} {definition}")

    @staticmethod
    def _input_fingerprint(smiles: str, notes: str | None, params_json: str | None) -> str:
        return fingerprint(
            settings={"notes": notes or "", "params_json": params_json or ""},
            inputs={"smiles": smiles},
            tool_version="input-v1",
        )

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

    def _read_ligand_rows(self, path: Path) -> list[dict[str, str | None]]:
        rows: list[dict[str, str | None]] = []
        seen: set[str] = set()
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or "id" not in reader.fieldnames or "smiles" not in reader.fieldnames:
                raise ValueError("Ligand CSV must contain stable 'id' and 'smiles' columns")
            for line_number, row in enumerate(reader, 2):
                parent_id = (row.get("id") or "").strip()
                smiles = (row.get("smiles") or "").strip()
                if not parent_id:
                    raise ValueError(f"Missing stable ligand id on CSV line {line_number}")
                if not smiles:
                    raise ValueError(f"Missing SMILES for ligand '{parent_id}' on CSV line {line_number}")
                if parent_id in seen:
                    raise ValueError(f"Duplicate ligand id '{parent_id}' in CSV")
                seen.add(parent_id)
                rows.append({
                    "parent_id": parent_id,
                    "source_id": parent_id,
                    "source_smiles": smiles,
                    "notes": (row.get("notes") or "").strip(),
                    "params_json": (row.get("params_json") or "").strip(),
                })
        if not rows:
            raise ValueError("Ligand CSV contains no usable rows")
        return rows

    def preview_ligand_sync(self, path: Path) -> LigandSyncResult:
        incoming = self._read_ligand_rows(path)
        with self.connection() as conn:
            existing = {row["parent_id"]: row for row in conn.execute("SELECT parent_id, active, input_fingerprint FROM parent_ligands").fetchall()}
        incoming_ids = {str(row["parent_id"]) for row in incoming}
        added = unchanged = changed = 0
        for row in incoming:
            parent_id = str(row["parent_id"])
            current = existing.get(parent_id)
            input_fp = self._input_fingerprint(str(row["source_smiles"]), row["notes"], row["params_json"])
            if current and int(current["active"]) == 1 and current["input_fingerprint"] == input_fp:
                unchanged += 1
            elif current:
                changed += 1
            else:
                added += 1
        archived = sum(1 for parent_id, row in existing.items() if int(row["active"]) == 1 and parent_id not in incoming_ids)
        return LigandSyncResult(added=added, unchanged=unchanged, changed=changed, archived=archived)

    def _invalidate_parent(self, conn: sqlite3.Connection, parent_id: str) -> None:
        conn.execute("UPDATE screening_results SET active=0 WHERE parent_id=?", (parent_id,))
        conn.execute("UPDATE molecular_states SET active=0, status='obsolete' WHERE parent_id=?", (parent_id,))
        conn.execute("""UPDATE conformers SET status='obsolete' WHERE state_id IN
            (SELECT state_id FROM molecular_states WHERE parent_id=?)""", (parent_id,))
        conn.execute("""UPDATE docking_runs SET is_current=0 WHERE state_id IN
            (SELECT state_id FROM molecular_states WHERE parent_id=?)""", (parent_id,))
        conn.execute("""UPDATE artifacts SET active=0 WHERE artifact_id IN
            (SELECT c.sdf_artifact_id FROM conformers c JOIN molecular_states s ON s.state_id=c.state_id WHERE s.parent_id=? AND c.sdf_artifact_id IS NOT NULL)
            OR artifact_id IN (SELECT d.log_artifact_id FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id WHERE s.parent_id=? AND d.log_artifact_id IS NOT NULL)
            OR artifact_id IN (SELECT d.raw_output_artifact_id FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id WHERE s.parent_id=? AND d.raw_output_artifact_id IS NOT NULL)""", (parent_id, parent_id, parent_id))

    def sync_ligands_csv(self, path: Path) -> LigandSyncResult:
        incoming = self._read_ligand_rows(path)
        incoming_ids = {str(row["parent_id"]) for row in incoming}
        with self.connection() as conn:
            existing = {row["parent_id"]: row for row in conn.execute("SELECT * FROM parent_ligands").fetchall()}
            counts = {"added": 0, "unchanged": 0, "changed": 0, "archived": 0}
            now = utc_now()
            for row in incoming:
                parent_id = str(row["parent_id"])
                input_fp = self._input_fingerprint(str(row["source_smiles"]), row["notes"], row["params_json"])
                current = existing.get(parent_id)
                if current and int(current["active"]) == 1 and current["input_fingerprint"] == input_fp:
                    counts["unchanged"] += 1
                    continue
                if current:
                    self._invalidate_parent(conn, parent_id)
                    conn.execute("""UPDATE parent_ligands SET source_id=?, source_smiles=?, notes=?, params_json=?,
                        canonical_source_smiles=NULL, parent_inchikey=NULL, parse_status=?, parse_reason=NULL,
                        active=1, input_fingerprint=?, archived_at=NULL WHERE parent_id=?""", (row["source_id"], row["source_smiles"], row["notes"], row["params_json"], StageStatus.PENDING, input_fp, parent_id))
                    counts["changed"] += 1
                else:
                    conn.execute("""INSERT INTO parent_ligands
                        (parent_id, source_id, source_smiles, notes, params_json, parse_status, created_at, active, input_fingerprint)
                        VALUES (?, ?, ?, ?, ?, ?, ?, 1, ?)""", (parent_id, row["source_id"], row["source_smiles"], row["notes"], row["params_json"], StageStatus.PENDING, now, input_fp))
                    counts["added"] += 1
            for parent_id, current in existing.items():
                if int(current["active"]) == 1 and parent_id not in incoming_ids:
                    self._invalidate_parent(conn, parent_id)
                    conn.execute("UPDATE parent_ligands SET active=0, archived_at=? WHERE parent_id=?", (now, parent_id))
                    counts["archived"] += 1
        return LigandSyncResult(**counts)

    def import_ligands_csv(self, path: Path) -> int:
        """Compatibility wrapper for callers that expect an import count."""
        result = self.sync_ligands_csv(path)
        return result.added + result.changed

    def clear_workflow_data(self) -> None:
        """Clear current ligand/stage data while retaining project settings."""
        with self.connection() as conn:
            conn.execute("PRAGMA foreign_keys = OFF")
            for table in ("docking_poses", "docking_runs", "conformers", "molecular_states", "screening_results", "stage_runs", "artifacts", "manual_reviews", "parent_ligands"):
                conn.execute(f"DELETE FROM {table}")
            conn.execute("PRAGMA foreign_keys = ON")

    def add_artifact(self, path: Path, artifact_type: str, stage: str, validation_status: str = "valid") -> str:
        path = path.resolve()
        relative = path.relative_to(self.root).as_posix()
        digest = file_sha256(path)
        with self.connection() as conn:
            existing = conn.execute("SELECT artifact_id FROM artifacts WHERE relative_path=?", (relative,)).fetchone()
            artifact_id = existing[0] if existing else str(uuid.uuid4())
            conn.execute("""INSERT INTO artifacts(artifact_id,relative_path,sha256,size_bytes,artifact_type,created_stage,validation_status,active,created_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, 1, ?)
                ON CONFLICT(relative_path) DO UPDATE SET sha256=excluded.sha256,size_bytes=excluded.size_bytes,
                artifact_type=excluded.artifact_type,created_stage=excluded.created_stage,validation_status=excluded.validation_status,active=1""", (artifact_id, relative, digest, path.stat().st_size, artifact_type, stage, validation_status, utc_now()))
            return artifact_id

    def recover_interrupted_runs(self) -> int:
        with self.connection() as conn:
            result = conn.execute("UPDATE stage_runs SET status=?, ended_at=?, reason=? WHERE status=?", (StageStatus.INTERRUPTED, utc_now(), "Application restarted before completion", StageStatus.RUNNING)).rowcount
            conn.execute("UPDATE docking_runs SET status=?, ended_at=? WHERE status=?", (StageStatus.INTERRUPTED, utc_now(), StageStatus.RUNNING))
            redocking = conn.execute("UPDATE redocking_runs SET status='interrupted', interrupted_at=?, error_message=? WHERE status='running'",
                                     (utc_now(), "Application restarted before completion")).rowcount
            conn.execute("UPDATE redocking_stages SET status='interrupted', finished_at=? WHERE status='running'", (utc_now(),))
            return result + redocking
