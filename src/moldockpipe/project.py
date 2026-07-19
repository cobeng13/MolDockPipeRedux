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


SCHEMA_VERSION = 6


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
        for relative in ("inputs", "inputs/receptors", "artifacts/sdf", "artifacts/pdbqt", "artifacts/docking", "logs", "exports", "For_PostDocking", "tools/vina"):
            (root / relative).mkdir(parents=True, exist_ok=True)
        config = {
            "name": name or root.name,
            "schema_version": SCHEMA_VERSION,
            "screening": {"policy": ScreeningPolicy.ANNOTATE_ONLY.value, "lipinski": True, "veber": True, "egan": False, "ghose": False, "boiled_egg": False},
            "molscrub": {"ph": 7.4, "enumerate_states": True, "max_states": 32, "stereochemistry_policy": "warn_and_continue", "fragment_policy": "manual_review"},
            "meeko": {"workers": 4},
            "postdock": {"mode": "split_and_sdf", "poses_per_compound": 3, "selected_parents": []},
            # A receptor is project-specific. New projects deliberately start
            # empty so a placeholder cannot accidentally participate in runs.
            "vina": {"profiles": []},
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
            CREATE TABLE IF NOT EXISTS workflow_runs (
                workflow_run_id TEXT PRIMARY KEY, workflow_type TEXT NOT NULL, status TEXT NOT NULL,
                settings_json TEXT NOT NULL DEFAULT '{}', environment_json TEXT NOT NULL DEFAULT '{}',
                started_at TEXT NOT NULL, finished_at TEXT, error_message TEXT
            );
            CREATE TABLE IF NOT EXISTS workflow_run_stages (
                workflow_run_id TEXT NOT NULL REFERENCES workflow_runs(workflow_run_id) ON DELETE CASCADE,
                stage_name TEXT NOT NULL, status TEXT NOT NULL, summary_json TEXT NOT NULL DEFAULT '{}',
                started_at TEXT, finished_at TEXT, error_message TEXT,
                PRIMARY KEY(workflow_run_id, stage_name)
            );
            CREATE TABLE IF NOT EXISTS provenance_events (
                event_id TEXT PRIMARY KEY, workflow_run_id TEXT REFERENCES workflow_runs(workflow_run_id) ON DELETE CASCADE,
                event_time TEXT NOT NULL, level TEXT NOT NULL, event_type TEXT NOT NULL, stage_name TEXT,
                entity_type TEXT, entity_id TEXT, receptor_profile_id TEXT, message TEXT,
                data_json TEXT NOT NULL DEFAULT '{}'
            );
            CREATE TABLE IF NOT EXISTS receptor_preparation_runs (
                preparation_run_id TEXT PRIMARY KEY, receptor_profile_id TEXT NOT NULL, receptor_name TEXT NOT NULL,
                status TEXT NOT NULL, source_path TEXT NOT NULL, source_sha256 TEXT NOT NULL,
                inventory_json TEXT NOT NULL DEFAULT '{}', decisions_json TEXT NOT NULL DEFAULT '{}',
                command_json TEXT NOT NULL DEFAULT '{}', artifacts_json TEXT NOT NULL DEFAULT '{}',
                warnings_json TEXT NOT NULL DEFAULT '[]', started_at TEXT NOT NULL, finished_at TEXT,
                error_message TEXT
            );
            CREATE TABLE IF NOT EXISTS report_snapshots (
                report_id TEXT PRIMARY KEY, title TEXT NOT NULL, scope_json TEXT NOT NULL,
                source_runs_json TEXT NOT NULL, generated_at TEXT NOT NULL,
                output_path TEXT, output_sha256 TEXT
            );
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
            CREATE TABLE IF NOT EXISTS redocking_queue_items (
                queue_id TEXT PRIMARY KEY, receptor_profile_id TEXT NOT NULL, receptor_profile_name TEXT NOT NULL,
                status TEXT NOT NULL, settings_json TEXT NOT NULL, redocking_run_id TEXT,
                current_stage TEXT, created_at TEXT NOT NULL, started_at TEXT, finished_at TEXT,
                error_message TEXT, cancel_requested INTEGER NOT NULL DEFAULT 0
            );
            CREATE INDEX IF NOT EXISTS idx_redocking_profile_status ON redocking_runs(receptor_profile_id, status);
            CREATE INDEX IF NOT EXISTS idx_redocking_queue_status_created ON redocking_queue_items(status, created_at);
            """)
            self._ensure_column(conn, "parent_ligands", "active", "INTEGER NOT NULL DEFAULT 1")
            self._ensure_column(conn, "parent_ligands", "input_fingerprint", "TEXT")
            self._ensure_column(conn, "parent_ligands", "archived_at", "TEXT")
            self._ensure_column(conn, "screening_results", "active", "INTEGER NOT NULL DEFAULT 1")
            self._ensure_column(conn, "molecular_states", "active", "INTEGER NOT NULL DEFAULT 1")
            self._ensure_column(conn, "docking_runs", "receptor_profile_id", "TEXT NOT NULL DEFAULT 'default'")
            self._ensure_column(conn, "docking_runs", "receptor_profile_name", "TEXT NOT NULL DEFAULT 'Default'")
            self._ensure_column(conn, "docking_runs", "reason", "TEXT")
            self._ensure_column(conn, "docking_runs", "workflow_run_id", "TEXT")
            conn.execute("""WITH ranked AS (
                    SELECT run_id, ROW_NUMBER() OVER (
                        PARTITION BY state_id, receptor_profile_id
                        ORDER BY CASE WHEN status='completed' THEN 0 ELSE 1 END,
                                 COALESCE(ended_at, started_at, '') DESC, run_id DESC) rank
                    FROM docking_runs WHERE is_current=1)
                UPDATE docking_runs SET is_current=0 WHERE run_id IN (SELECT run_id FROM ranked WHERE rank>1)""")
            conn.execute("CREATE INDEX IF NOT EXISTS idx_docking_runs_profile_state ON docking_runs(receptor_profile_id, state_id, is_current)")
            conn.execute("CREATE INDEX IF NOT EXISTS idx_docking_runs_workflow ON docking_runs(workflow_run_id)")
            conn.execute("CREATE INDEX IF NOT EXISTS idx_workflow_runs_type_status ON workflow_runs(workflow_type, status)")
            conn.execute("CREATE INDEX IF NOT EXISTS idx_receptor_preparation_profile ON receptor_preparation_runs(receptor_profile_id, started_at)")
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

    def remove_receptor_profile(self, profile_id: str) -> bool:
        """Remove one receptor profile and only its receptor-owned files.

        Docking records are deliberately retained in the database for audit
        purposes, but become inactive because their profile is no longer in
        project.yml. Ligand inputs and shared artifacts are never touched.
        """
        profile_id = str(profile_id)
        profiles = self.get_receptor_profiles(include_archived=True)
        retained = [profile for profile in profiles if str(profile.get("id")) != profile_id]
        if len(retained) == len(profiles):
            return False
        with self.connection() as conn:
            active = conn.execute("""SELECT EXISTS(SELECT 1 FROM redocking_queue_items
                WHERE receptor_profile_id=? AND status='running')""", (profile_id,)).fetchone()[0]
            if active:
                raise RuntimeError("Cannot delete a receptor while its validation is running")
            conn.execute("""UPDATE redocking_queue_items SET status='cancelled',finished_at=?,
                error_message='Receptor profile was deleted',cancel_requested=1
                WHERE receptor_profile_id=? AND status='queued'""", (utc_now(), profile_id))
        self.save_receptor_profiles(retained)
        parent = (self.root / "inputs" / "receptors").resolve()
        target = (parent / profile_id).resolve()
        if target.parent != parent:
            raise ValueError("Invalid receptor profile identifier")
        if target.exists():
            shutil.rmtree(target)
        return True

    def enqueue_redocking(self, profile: dict[str, Any], settings: dict[str, Any]) -> str:
        """Persist one receptor-validation request for FIFO background execution."""
        profile_id = str(profile["id"])
        with self.connection() as conn:
            duplicate = conn.execute("""SELECT queue_id FROM redocking_queue_items
                WHERE receptor_profile_id=? AND status IN ('queued','running') LIMIT 1""", (profile_id,)).fetchone()
            if duplicate:
                raise ValueError("This receptor already has a queued or running validation.")
            queue_id = uuid.uuid4().hex
            conn.execute("""INSERT INTO redocking_queue_items
                (queue_id,receptor_profile_id,receptor_profile_name,status,settings_json,created_at)
                VALUES (?, ?, ?, 'queued', ?, ?)""",
                (queue_id, profile_id, str(profile.get("name", profile_id)),
                 json.dumps(settings, sort_keys=True), utc_now()))
        return queue_id

    def list_redocking_queue(self) -> list[dict[str, Any]]:
        with self.connection() as conn:
            rows = conn.execute("""SELECT * FROM redocking_queue_items
                ORDER BY CASE status WHEN 'running' THEN 0 WHEN 'queued' THEN 1 ELSE 2 END,
                         CASE WHEN status IN ('running','queued') THEN created_at END ASC,
                         CASE WHEN status NOT IN ('running','queued') THEN COALESCE(finished_at,created_at) END DESC""").fetchall()
        result = []
        for row in rows:
            item = dict(row)
            item["settings"] = json.loads(item.pop("settings_json"))
            result.append(item)
        return result

    def claim_next_redocking(self) -> dict[str, Any] | None:
        """Atomically claim the oldest queued validation."""
        with self.connection() as conn:
            row = conn.execute("""SELECT * FROM redocking_queue_items WHERE status='queued'
                ORDER BY created_at,rowid LIMIT 1""").fetchone()
            if not row:
                return None
            run_id = str(row["redocking_run_id"] or uuid.uuid4().hex)
            now = utc_now()
            updated = conn.execute("""UPDATE redocking_queue_items SET status='running',
                redocking_run_id=?,started_at=COALESCE(started_at,?),finished_at=NULL,
                error_message=NULL,cancel_requested=0 WHERE queue_id=? AND status='queued'""",
                (run_id, now, row["queue_id"])).rowcount
            if not updated:
                return None
            item = dict(row)
            item.update(status="running", redocking_run_id=run_id, started_at=item.get("started_at") or now)
            item["settings"] = json.loads(item.pop("settings_json"))
            return item

    def update_redocking_queue_stage(self, queue_id: str, stage: str) -> None:
        with self.connection() as conn:
            conn.execute("UPDATE redocking_queue_items SET current_stage=? WHERE queue_id=? AND status='running'",
                         (stage, queue_id))

    def finish_redocking_queue_item(self, queue_id: str, status: str, error: str | None = None) -> None:
        if status not in {"completed", "failed", "cancelled"}:
            raise ValueError(f"Invalid redocking queue terminal status: {status}")
        with self.connection() as conn:
            conn.execute("""UPDATE redocking_queue_items SET status=?,current_stage=NULL,finished_at=?,
                error_message=? WHERE queue_id=?""", (status, utc_now(), error, queue_id))

    def cancel_redocking_queue_item(self, queue_id: str) -> bool:
        with self.connection() as conn:
            row = conn.execute("SELECT status FROM redocking_queue_items WHERE queue_id=?", (queue_id,)).fetchone()
            if not row or row["status"] not in {"queued", "running"}:
                return False
            if row["status"] == "queued":
                conn.execute("""UPDATE redocking_queue_items SET status='cancelled',cancel_requested=1,
                    finished_at=?,error_message='Cancelled by user' WHERE queue_id=?""", (utc_now(), queue_id))
            else:
                conn.execute("UPDATE redocking_queue_items SET cancel_requested=1 WHERE queue_id=?", (queue_id,))
            return True

    def retry_redocking_queue_item(self, queue_id: str) -> bool:
        with self.connection() as conn:
            return bool(conn.execute("""UPDATE redocking_queue_items SET status='queued',current_stage=NULL,
                finished_at=NULL,error_message=NULL,cancel_requested=0 WHERE queue_id=?
                AND status IN ('failed','cancelled')""", (queue_id,)).rowcount)

    def requeue_running_redocking(self, queue_id: str, reason: str) -> bool:
        with self.connection() as conn:
            return bool(conn.execute("""UPDATE redocking_queue_items SET status='queued',current_stage=NULL,
                error_message=?,cancel_requested=0 WHERE queue_id=? AND status='running'""",
                (reason, queue_id)).rowcount)

    def clear_finished_redocking_queue(self) -> int:
        with self.connection() as conn:
            return conn.execute("DELETE FROM redocking_queue_items WHERE status IN ('completed','cancelled')").rowcount

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

    def _environment_snapshot(self, conn: sqlite3.Connection) -> dict[str, Any]:
        row = conn.execute("SELECT value_json FROM project_settings WHERE key='environment'").fetchone()
        environment = json.loads(row[0]) if row else {
            "os": platform.platform(), "architecture": platform.machine(), "python": sys.version,
        }
        environment["tools"] = [dict(tool) for tool in conn.execute(
            "SELECT tool_name,version,location,sha256,recorded_at FROM tool_installations ORDER BY tool_name"
        ).fetchall()]
        return environment

    def start_workflow_run(self, workflow_type: str, settings: dict[str, Any]) -> str:
        workflow_run_id = uuid.uuid4().hex
        with self.connection() as conn:
            environment = self._environment_snapshot(conn)
            conn.execute("""INSERT INTO workflow_runs
                (workflow_run_id,workflow_type,status,settings_json,environment_json,started_at)
                VALUES (?, ?, ?, ?, ?, ?)""", (workflow_run_id, workflow_type, StageStatus.RUNNING.value,
                json.dumps(settings, sort_keys=True), json.dumps(environment, sort_keys=True), utc_now()))
        return workflow_run_id

    def update_workflow_stage(self, workflow_run_id: str, stage_name: str, status: str,
                              *, summary: dict[str, Any] | None = None, error: str | None = None) -> None:
        now = utc_now()
        started_at = now if status == StageStatus.RUNNING.value else None
        finished_at = now if status in {StageStatus.COMPLETED.value, StageStatus.FAILED.value,
                                        StageStatus.CANCELLED.value, StageStatus.INTERRUPTED.value} else None
        with self.connection() as conn:
            conn.execute("""INSERT INTO workflow_run_stages
                (workflow_run_id,stage_name,status,summary_json,started_at,finished_at,error_message)
                VALUES (?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(workflow_run_id,stage_name) DO UPDATE SET
                status=excluded.status,summary_json=excluded.summary_json,
                started_at=COALESCE(workflow_run_stages.started_at,excluded.started_at),
                finished_at=excluded.finished_at,error_message=excluded.error_message""",
                (workflow_run_id, stage_name, status, json.dumps(summary or {}, sort_keys=True),
                 started_at, finished_at, error))

    def finish_workflow_run(self, workflow_run_id: str, status: str,
                            *, error: str | None = None) -> None:
        with self.connection() as conn:
            conn.execute("UPDATE workflow_runs SET status=?,finished_at=?,error_message=? WHERE workflow_run_id=?",
                         (status, utc_now(), error, workflow_run_id))

    def record_provenance_event(self, *, event_type: str, level: str = "INFO",
                                workflow_run_id: str | None = None, stage_name: str | None = None,
                                entity_type: str | None = None, entity_id: str | None = None,
                                receptor_profile_id: str | None = None, message: str | None = None,
                                data: dict[str, Any] | None = None) -> str:
        event_id = uuid.uuid4().hex
        with self.connection() as conn:
            conn.execute("INSERT INTO provenance_events VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                         (event_id, workflow_run_id, utc_now(), level, event_type, stage_name,
                          entity_type, entity_id, receptor_profile_id, message,
                          json.dumps(data or {}, sort_keys=True)))
        return event_id

    def start_receptor_preparation(self, *, profile_id: str, receptor_name: str,
                                   source_path: Path, source_sha256: str,
                                   inventory: dict[str, Any], decisions: dict[str, Any]) -> str:
        preparation_run_id = uuid.uuid4().hex
        with self.connection() as conn:
            conn.execute("""INSERT INTO receptor_preparation_runs
                (preparation_run_id,receptor_profile_id,receptor_name,status,source_path,source_sha256,
                 inventory_json,decisions_json,started_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (preparation_run_id, profile_id, receptor_name, StageStatus.RUNNING.value,
                 str(source_path), source_sha256, json.dumps(inventory, sort_keys=True),
                 json.dumps(decisions, sort_keys=True), utc_now()))
        return preparation_run_id

    def finish_receptor_preparation(self, preparation_run_id: str, status: str, *,
                                    command: dict[str, Any] | None = None,
                                    artifacts: dict[str, Any] | None = None,
                                    warnings: list[str] | None = None,
                                    error: str | None = None) -> None:
        with self.connection() as conn:
            conn.execute("""UPDATE receptor_preparation_runs SET status=?,command_json=?,artifacts_json=?,
                warnings_json=?,finished_at=?,error_message=? WHERE preparation_run_id=?""",
                (status, json.dumps(command or {}, sort_keys=True), json.dumps(artifacts or {}, sort_keys=True),
                 json.dumps(warnings or []), utc_now(), error, preparation_run_id))

    def record_report_snapshot(self, *, report_id: str, title: str, scope: dict[str, Any],
                               source_runs: dict[str, Any], output_path: Path,
                               output_sha256: str) -> None:
        relative = output_path.resolve().relative_to(self.root).as_posix()
        with self.connection() as conn:
            conn.execute("""INSERT INTO report_snapshots
                (report_id,title,scope_json,source_runs_json,generated_at,output_path,output_sha256)
                VALUES (?, ?, ?, ?, ?, ?, ?)""",
                (report_id, title, json.dumps(scope, sort_keys=True),
                 json.dumps(source_runs, sort_keys=True), utc_now(), relative, output_sha256))

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
            for table in ("provenance_events", "workflow_run_stages", "workflow_runs", "report_snapshots",
                          "docking_poses", "docking_runs", "conformers", "molecular_states",
                          "screening_results", "stage_runs", "artifacts", "manual_reviews", "parent_ligands"):
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
            workflows = conn.execute("UPDATE workflow_runs SET status='interrupted',finished_at=?,error_message=? WHERE status='running'",
                                     (utc_now(), "Application restarted before completion")).rowcount
            conn.execute("UPDATE workflow_run_stages SET status='interrupted',finished_at=?,error_message=? WHERE status='running'",
                         (utc_now(), "Application restarted before completion"))
            queued = conn.execute("""UPDATE redocking_queue_items SET status='queued',current_stage=NULL,
                error_message='Application restarted; validation will resume from reusable stages',cancel_requested=0
                WHERE status='running'""").rowcount
            return result + redocking + workflows + queued
