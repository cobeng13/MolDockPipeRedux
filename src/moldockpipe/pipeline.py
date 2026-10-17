from __future__ import annotations

import json
import os
import uuid
from collections import defaultdict
from pathlib import Path
from typing import Collection, Callable

from .fingerprints import file_sha256, fingerprint
from .models import ProgressCallback, ProgressEvent, ScreeningPolicy, StageStatus
from .project import ProjectRepository, utc_now
from .parallel import bounded_futures, logical_cpu_count, worker_count
from .execution_lock import exclusive_project
from .services.molscrub import MolScrubService
from .services.meeko import MeekoService
from .services.screening import screen_smiles
from .services.validation import validate_pdbqt
from .services.vina import VinaDockingBackend, find_vina_executable, find_vina_split_executable
from .services.vinagpu import VinaGpuDockingBackend, find_vinagpu_executable, vinagpu_tool_hash
from .services.postdock import PostDockService, receptor_export_filename


PLANNING_BATCH_SIZE = 256


def _resolve_destination(path: Path, cache: dict[Path, Path]) -> Path:
    """Resolve existing paths normally; share ancestor work for missing outputs.

    A missing output need not repeatedly resolve the same missing directories.
    Existing entries (including symlinks, junctions and Windows aliases) still
    use Path.resolve(), preserving physical destination collision checks.
    The cache belongs only to the read-only preflight, before any worker writes.
    """
    if os.name == "nt":
        # Windows realpath normalizes '..' before following junctions/links;
        # POSIX realpath follows each component before interpreting '..'.
        path = Path(os.path.normpath(path))
    if path in cache:
        return cache[path]
    try:
        path.lstat()
    except FileNotFoundError:
        if path.parent == path:
            resolved = path.resolve()
        else:
            parent = _resolve_destination(path.parent, cache)
            resolved = parent.parent if path.name == ".." else parent / path.name
    except OSError:
        # Retain pathlib's platform-specific error handling for other failures.
        resolved = path.resolve()
    else:
        resolved = path.resolve()
    cache[path] = resolved
    return resolved


def _docking_parameters(settings: dict[str, object], backend: str) -> dict[str, object]:
    """Fingerprint only arguments used by the selected docking executable."""
    defaults = {"num_modes": 9, "energy_range": 3, "seed": 0}
    if backend == "vinagpu":
        defaults.update(thread=8000, search_depth=32)
    else:
        defaults.update(exhaustiveness=8, cpu=0)
    names = ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z", *defaults)
    # JSON distinguishes 20 from 20.0, although Vina treats them identically.
    return {"backend": backend, **{name: float(settings.get(name, defaults.get(name)))
            for name in names if settings.get(name, defaults.get(name)) is not None}}


class PipelineRunner:
    """Application service for small, independently resumable stage operations."""

    def __init__(self, repository: ProjectRepository, progress: ProgressCallback | None = None,
                 cancelled: Callable[[], bool] | None = None):
        self.repository = repository
        self.progress = progress
        self.cancelled = cancelled

    def _check_cancelled(self) -> None:
        if self.cancelled and self.cancelled():
            raise InterruptedError("Pipeline stopped safely")

    def _emit(self, event: str, stage: str, item_id: str | None = None, index: int = 0, total: int = 0, **kwargs: object) -> None:
        if self.progress:
            self.progress(ProgressEvent(event, stage, item_id, index, total, **kwargs))

    def _artifact_matches(self, relative_path: str | None, sha256: str | None) -> bool:
        if not relative_path or not sha256:
            return False
        path = self.repository.root / relative_path
        return path.is_file() and file_sha256(path) == sha256

    def _molscrub_plan_rows(self, rows):
        """Read history/states in bounded batches without changing reuse rules."""
        for offset in range(0, len(rows), PLANNING_BATCH_SIZE):
            self._check_cancelled()
            batch = rows[offset:offset + PLANNING_BATCH_SIZE]
            ids = [row["parent_id"] for row in batch]
            placeholders = ",".join("?" for _ in ids)
            with self.repository.connection() as conn:
                latest = {row["entity_id"]: row for row in conn.execute(f"""
                    SELECT entity_id,fingerprint,status,reason FROM stage_runs WHERE rowid IN (
                        SELECT MAX(rowid) FROM stage_runs WHERE stage_name='molscrub'
                        AND entity_type='parent' AND entity_id IN ({placeholders}) GROUP BY entity_id)
                    """, ids)}
                states = defaultdict(list)
                for state in conn.execute(f"""SELECT s.parent_id,s.state_structure_hash,s.generation_fingerprint,
                    a.relative_path,a.sha256 FROM molecular_states s LEFT JOIN conformers c ON c.state_id=s.state_id
                    LEFT JOIN artifacts a ON a.artifact_id=c.sdf_artifact_id
                    WHERE s.parent_id IN ({placeholders}) AND s.active=1""", ids):
                    states[state["parent_id"]].append(state)
            for row in batch:
                yield row, latest.get(row["parent_id"]), states[row["parent_id"]]

    def _meeko_plan_rows(self, rows):
        for offset in range(0, len(rows), PLANNING_BATCH_SIZE):
            self._check_cancelled()
            batch = rows[offset:offset + PLANNING_BATCH_SIZE]
            paths = [(self.repository.root / "artifacts" / "pdbqt" / row["parent_id"] /
                      row["state_id"] / "ligand.pdbqt").relative_to(self.repository.root).as_posix()
                     for row in batch]
            with self.repository.connection() as conn:
                artifacts = {artifact["relative_path"]: artifact for artifact in conn.execute(
                    f"SELECT relative_path,sha256 FROM artifacts WHERE active=1 AND relative_path IN ({','.join('?' for _ in paths)})",
                    paths)}
            for row, path in zip(batch, paths):
                yield row, artifacts.get(path)

    def _docking_plan_rows(self, states, profile_id):
        seen = set()
        for offset in range(0, len(states), PLANNING_BATCH_SIZE):
            self._check_cancelled()
            batch = states[offset:offset + PLANNING_BATCH_SIZE]
            ids = [state["state_id"] for state in batch]
            paths = [(self.repository.root / "artifacts" / "pdbqt" / state["parent_id"] /
                      state["state_id"] / "ligand.pdbqt").relative_to(self.repository.root).as_posix()
                     for state in batch]
            placeholders = ",".join("?" for _ in ids)
            with self.repository.connection() as conn:
                prepared = {row["relative_path"]: row for row in conn.execute(
                    f"SELECT relative_path,sha256 FROM artifacts WHERE active=1 AND relative_path IN ({placeholders})", paths)}
                candidates = defaultdict(list)
                for run in conn.execute(f"""SELECT d.*,raw.relative_path AS raw_path,raw.sha256 AS raw_sha256
                    FROM docking_runs d JOIN artifacts raw ON raw.artifact_id=d.raw_output_artifact_id
                    WHERE d.state_id IN ({placeholders}) AND d.receptor_profile_id=? AND d.status='completed'
                    AND EXISTS (SELECT 1 FROM docking_poses p WHERE p.run_id=d.run_id)
                    ORDER BY d.is_current DESC,d.rowid DESC""", [*ids, profile_id]):
                    candidates[run["state_id"]].append(run)
            for state, path in zip(batch, paths):
                # Old/imported projects may contain several conformers per
                # state. A preceding occurrence may just have written a run;
                # retain the original live lookup for duplicate occurrences.
                if state["state_id"] in seen:
                    yield state, None, None
                else:
                    seen.add(state["state_id"])
                    yield state, prepared.get(path), candidates[state["state_id"]]

    def _reusable_docking(self, state_id: str, profile_id: str, settings_fp: str,
                          parameters: dict[str, object], receptor_hash: str,
                          tool_hash: str, backend: str, ligand: Path, ligand_hash: str,
                          *, candidates=None, prepared=None):
        if candidates is None:
            with self.repository.connection() as conn:
                prepared = conn.execute("SELECT sha256 FROM artifacts WHERE relative_path=? AND active=1",
                                       (ligand.relative_to(self.repository.root).as_posix(),)).fetchone()
                candidates = conn.execute("""SELECT d.*, raw.relative_path AS raw_path,
                    raw.sha256 AS raw_sha256 FROM docking_runs d
                    JOIN artifacts raw ON raw.artifact_id=d.raw_output_artifact_id
                    WHERE d.state_id=? AND d.receptor_profile_id=? AND d.status='completed'
                    AND EXISTS (SELECT 1 FROM docking_poses p WHERE p.run_id=d.run_id)
                    ORDER BY d.is_current DESC, d.rowid DESC""", (state_id, profile_id)).fetchall()
        for run in candidates:
            matches = run["settings_fingerprint"] == settings_fp
            if (not matches and run["is_current"] and run["receptor_hash"] == receptor_hash
                    and prepared and prepared["sha256"] == ligand_hash):
                # Keep compatible pre-fingerprint-upgrade results reusable. Old
                # non-current runs lack the ligand hash needed to prove that
                # regenerated coordinates are the same input, so cannot qualify.
                command = json.loads(run["command_json"] or "{}")
                previous = command.get("settings", {})
                matches = ("inputs" not in command
                    and _docking_parameters(previous, backend) == parameters
                    and fingerprint(settings=previous, inputs={"receptor": receptor_hash},
                                    tool_version=tool_hash) == run["settings_fingerprint"])
            if matches and self._artifact_matches(run["raw_path"], run["raw_sha256"]):
                return run
        return None

    @exclusive_project
    def run_screening(self) -> tuple[int, int]:
        self._check_cancelled()
        settings = self.repository.get_settings().get("screening", {})
        policy = ScreeningPolicy(settings.get("policy", ScreeningPolicy.ANNOTATE_ONLY.value))
        enabled = {name: bool(settings.get(name, False)) for name in ("lipinski", "veber", "egan", "ghose", "boiled_egg")}
        try:
            import rdkit
            rdkit_version = rdkit.__version__
        except ImportError as exc:
            raise RuntimeError("RDKit must be installed before screening") from exc
        completed = failed = 0
        with self.repository.connection() as conn:
            parents = conn.execute("""SELECT p.parent_id,p.source_smiles,
                s.fingerprint AS previous_fingerprint,s.status AS previous_status
                FROM parent_ligands p LEFT JOIN screening_results s ON s.parent_id=p.parent_id AND s.active=1
                ORDER BY p.parent_id""").fetchall()
        self._emit("stage_started", "screening", total=len(parents))
        skipped = 0
        for index, parent in enumerate(parents, 1):
            self._check_cancelled()
            self._emit("item_started", "screening", parent["parent_id"], index, len(parents))
            run_fingerprint = fingerprint(settings={"enabled": enabled, "policy": policy.value}, inputs={"smiles": parent["source_smiles"]}, tool_version=rdkit_version)
            if parent["previous_fingerprint"] == run_fingerprint and parent["previous_status"] in (StageStatus.COMPLETED, StageStatus.FAILED):
                skipped += 1
                if parent["previous_status"] == StageStatus.COMPLETED:
                    completed += 1
                else:
                    failed += 1
                self._emit("item_skipped", "screening", parent["parent_id"], index, len(parents), skipped=skipped, message="Matching screening result reused")
                continue
            result = screen_smiles(parent["source_smiles"], enabled, policy)
            run_id = str(uuid.uuid4())
            with self.repository.connection() as conn:
                conn.execute("INSERT INTO stage_runs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", (run_id, "screening", "parent", parent["parent_id"], run_fingerprint, StageStatus.RUNNING, utc_now(), None, None))
                conn.execute("UPDATE parent_ligands SET canonical_source_smiles=?, parent_inchikey=?, parse_status=?, parse_reason=? WHERE parent_id=?", (result.canonical_smiles, result.inchikey, StageStatus.COMPLETED if result.decision != "fail" else StageStatus.FAILED, result.reason, parent["parent_id"]))
                conn.execute("""INSERT OR REPLACE INTO screening_results
                    (parent_id, fingerprint, status, descriptors_json, rules_json, decision, reason, created_at, active)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, 1)
                    """, (parent["parent_id"], run_fingerprint, StageStatus.COMPLETED if result.decision != "fail" else StageStatus.FAILED, json.dumps(result.descriptors), json.dumps(result.rules), result.decision, result.reason, utc_now()))
                conn.execute("UPDATE stage_runs SET status=?, ended_at=?, reason=? WHERE stage_run_id=?", (StageStatus.COMPLETED if result.decision != "fail" else StageStatus.FAILED, utc_now(), result.reason, run_id))
            if result.decision == "fail":
                failed += 1
                self._emit("item_failed", "screening", parent["parent_id"], index, len(parents), failed=failed, message=result.reason)
            else:
                completed += 1
                self._emit("item_succeeded", "screening", parent["parent_id"], index, len(parents), succeeded=completed)
        self.repository.record_tool("rdkit", rdkit_version, None)
        self._emit("stage_completed", "screening", total=len(parents), succeeded=completed, skipped=skipped, failed=failed)
        return completed, failed

    @exclusive_project
    def run_molscrub(self, *, retry_failed: bool = False) -> tuple[int, int]:
        self._check_cancelled()
        project = self.repository.get_settings()
        settings = project.get("molscrub", {})
        screening = project.get("screening", {})
        policy = screening.get("policy", ScreeningPolicy.ANNOTATE_ONLY.value)
        with self.repository.connection() as conn:
            rows = conn.execute("""SELECT p.parent_id, p.source_smiles, s.decision FROM parent_ligands p
                LEFT JOIN screening_results s ON s.parent_id=p.parent_id AND s.active=1
                ORDER BY p.parent_id""").fetchall()
        allowed = [row for row in rows if row["decision"] not in ("fail", None) or policy == ScreeningPolicy.ANNOTATE_ONLY.value]
        service = MolScrubService()
        version = getattr(self, "_molscrub_tool_version", None)
        if version is None:
            self._emit("stage_preparing", "molscrub", message="Loading MolScrub dependencies (30-second limit; Stop is available)")
            version = service.version(cancelled=self.cancelled)
        self._molscrub_tool_version = version
        created = failed = skipped = 0
        completed = 0
        workers = worker_count(settings, len(allowed))
        # Worker count affects scheduling, not chemical identity or failure
        # eligibility. Retain compatibility with fingerprints from serial runs.
        cache_settings = {key: value for key, value in settings.items() if key != "workers"}
        generation_settings = {key: value for key, value in cache_settings.items() if key != "timeout_seconds"}
        self._emit("stage_started", "molscrub", total=len(allowed), message=f"Using {workers} parallel workers")

        def jobs():
            nonlocal created, failed, skipped, completed
            for row, previous, existing in self._molscrub_plan_rows(allowed):
                if self.cancelled and self.cancelled():
                    return
                run_fingerprint = fingerprint(settings=cache_settings,
                    inputs={"smiles": row["source_smiles"]}, tool_version=version)
                if not retry_failed and previous and previous["status"] == StageStatus.FAILED and previous["fingerprint"] == run_fingerprint:
                    failed += 1
                    skipped += 1
                    completed += 1
                    self._emit("item_skipped", "molscrub", row["parent_id"], completed, len(allowed),
                               failed=failed, skipped=skipped, message=f"Known failure skipped: {previous['reason']}")
                    continue
                if (previous and previous["status"] == StageStatus.COMPLETED and existing and all(
                    state["generation_fingerprint"] in {
                        fingerprint(settings=candidate, inputs={"state": state["state_structure_hash"]}, tool_version=version)
                        for candidate in (generation_settings, cache_settings)}
                    and self._artifact_matches(state["relative_path"], state["sha256"])
                    for state in existing
                )):
                    created += len(existing)
                    skipped += 1
                    completed += 1
                    self._emit("item_skipped", "molscrub", row["parent_id"], completed, len(allowed),
                               skipped=skipped, message=f"{len(existing)} existing states reused")
                    continue
                run_id = str(uuid.uuid4())
                with self.repository.connection() as conn:
                    conn.execute("INSERT INTO stage_runs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                                 (run_id, "molscrub", "parent", row["parent_id"], run_fingerprint, StageStatus.RUNNING, utc_now(), None, None))
                self._emit("item_started", "molscrub", row["parent_id"], completed, len(allowed))
                yield row, run_id

        def generate(job):
            row, _ = job
            return service.generate(row["source_smiles"], ph=float(settings.get("ph", 7.4)),
                enumerate_states=bool(settings.get("enumerate_states", True)),
                max_states=int(settings.get("max_states", 32)),
                timeout_seconds=float(settings.get("timeout_seconds", 15)))

        def discard(job):
            _, run_id = job
            with self.repository.connection() as conn:
                conn.execute("UPDATE stage_runs SET status=?, ended_at=?, reason=? WHERE stage_run_id=?",
                             (StageStatus.CANCELLED, utc_now(), "Stopped before generation started", run_id))

        for (row, run_id), future in bounded_futures(generate, jobs(), workers,
                cancelled=self.cancelled, name="molscrub", on_discard=discard):
            completed += 1
            try:
                states, truncated = future.result()
                if not states:
                    raise ValueError("MolScrub returned no molecular states")
                self._persist_states(row["parent_id"], states, truncated, settings)
                with self.repository.connection() as conn:
                    conn.execute("UPDATE stage_runs SET status=?, ended_at=? WHERE stage_run_id=?", (StageStatus.COMPLETED, utc_now(), run_id))
                    conn.execute("UPDATE parent_ligands SET parse_reason=NULL WHERE parent_id=? AND parse_reason LIKE 'MolScrub failed:%'", (row["parent_id"],))
                created += len(states)
                self._emit("item_succeeded", "molscrub", row["parent_id"], completed, len(allowed), succeeded=created)
            except Exception as exc:
                failed += 1
                with self.repository.connection() as conn:
                    conn.execute("UPDATE parent_ligands SET parse_reason=? WHERE parent_id=?", (f"MolScrub failed: {exc}", row["parent_id"]))
                    conn.execute("UPDATE stage_runs SET status=?, ended_at=?, reason=? WHERE stage_run_id=?", (StageStatus.FAILED, utc_now(), str(exc), run_id))
                    # Partial or stale states from a failed parent must not enter downstream stages.
                    conn.execute("UPDATE molecular_states SET active=0 WHERE parent_id=?", (row["parent_id"],))
                    conn.execute("UPDATE docking_runs SET is_current=0 WHERE state_id IN (SELECT state_id FROM molecular_states WHERE parent_id=?)", (row["parent_id"],))
                self._emit("item_failed", "molscrub", row["parent_id"], completed, len(allowed), failed=failed, message=str(exc))
        self._emit("stage_completed", "molscrub", total=len(allowed), succeeded=created, skipped=skipped, failed=failed)
        return created, failed

    def _persist_states(self, parent_id: str, states: list[object], truncated: bool, settings: dict[str, object]) -> None:
        from rdkit import Chem
        version = getattr(self, "_molscrub_tool_version", None)
        if version is None:
            version = MolScrubService.version(cancelled=self.cancelled)
        with self.repository.connection() as conn:
            conn.execute("UPDATE molecular_states SET active=0 WHERE parent_id=?", (parent_id,))
            conn.execute("""UPDATE docking_runs SET is_current=0 WHERE state_id IN
                (SELECT state_id FROM molecular_states WHERE parent_id=?)""", (parent_id,))
        for state in states:
            state_id = f"{parent_id}-{state.structure_hash[:16]}"
            generation_settings = {key: value for key, value in settings.items() if key not in {"timeout_seconds", "workers"}}
            state_fingerprint = fingerprint(settings=generation_settings, inputs={"state": state.structure_hash}, tool_version=version)
            state_dir = self.repository.root / "artifacts" / "sdf" / parent_id / state_id
            state_dir.mkdir(parents=True, exist_ok=True)
            sdf_path = state_dir / "input_conformer_0001.sdf"
            writer = Chem.SDWriter(str(sdf_path))
            writer.write(state.molecule)
            writer.close()
            with self.repository.connection() as conn:
                conn.execute("""INSERT INTO molecular_states
                    (state_id,parent_id,state_smiles,state_isomeric_smiles,state_inchikey,formal_charge,tautomer_index,protomer_index,state_structure_hash,generation_fingerprint,status,reason,enumeration_truncated,created_at)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ON CONFLICT(state_id) DO UPDATE SET generation_fingerprint=excluded.generation_fingerprint,
                    status=excluded.status,reason=NULL,enumeration_truncated=excluded.enumeration_truncated,
                    created_at=excluded.created_at,active=1""", (state_id, parent_id, state.state_smiles, state.state_isomeric_smiles, state.state_inchikey, state.formal_charge, state.tautomer_index, state.protomer_index, state.structure_hash, state_fingerprint, "sdf_ready", None, int(truncated), utc_now()))
            artifact_id = self.repository.add_artifact(sdf_path, "state_sdf", "molscrub")
            with self.repository.connection() as conn:
                conn.execute("INSERT OR REPLACE INTO conformers VALUES (?, ?, ?, ?, ?, ?)", (f"{state_id}-conf-1", state_id, 1, artifact_id, "ready", utc_now()))

    @exclusive_project
    def run_meeko(self) -> tuple[int, int]:
        """Prepare every retained SDF state into a state-specific ligand PDBQT."""
        self._check_cancelled()
        prepared = failed = 0
        settings = self.repository.get_settings().get("meeko", {})
        with self.repository.connection() as conn:
            all_states = conn.execute("SELECT COUNT(*) FROM molecular_states WHERE active=1").fetchone()[0]
            rows = conn.execute("""SELECT s.state_id, s.parent_id, c.conformer_id, c.status AS conformer_status,a.relative_path
                FROM molecular_states s JOIN conformers c ON c.state_id=s.state_id
                JOIN artifacts a ON a.artifact_id=c.sdf_artifact_id WHERE s.active=1 AND a.active=1
                ORDER BY s.parent_id, s.state_id""").fetchall()
        self._emit("stage_started", "meeko", total=len(rows), message="Checking existing ligand outputs")
        reusable = set()
        destinations: dict[Path, str] = {}
        resolved_paths: dict[Path, Path] = {}
        for row, artifact in self._meeko_plan_rows(rows):
            self._check_cancelled()
            pdbqt = self.repository.root / "artifacts" / "pdbqt" / row["parent_id"] / row["state_id"] / "ligand.pdbqt"
            log = self.repository.root / "logs" / "meeko" / f"{row['state_id']}.log"
            for destination in (_resolve_destination(pdbqt, resolved_paths), _resolve_destination(log, resolved_paths)):
                if destination in destinations:
                    raise ValueError("Multiple Meeko preparation jobs target the same output file: "
                                     f"{destinations[destination]} and {row['conformer_id']}. "
                                     "Each state must have one active conformer.")
                destinations[destination] = row["conformer_id"]
            if row["conformer_status"] == "pdbqt_ready" and artifact and self._artifact_matches(artifact["relative_path"], artifact["sha256"]):
                reusable.add(row["conformer_id"])
        rows = [row for row in rows if row["conformer_id"] not in reusable]
        del resolved_paths, destinations
        workers = worker_count(settings, len(rows))
        self._emit("planning_completed", "meeko", total=len(rows), skipped=max(0, all_states - len(rows)),
                   message=f"Using {workers} parallel workers; {len(rows)} ligand preparations pending")
        skipped = max(0, all_states - len(rows))
        def prepare(row: object) -> tuple[object, Path, Path, str | None]:
            sdf = self.repository.root / row["relative_path"]
            output_dir = self.repository.root / "artifacts" / "pdbqt" / row["parent_id"] / row["state_id"]
            pdbqt = output_dir / "ligand.pdbqt"
            log = self.repository.root / "logs" / "meeko" / f"{row['state_id']}.log"
            try:
                MeekoService().prepare_ligand(sdf, pdbqt, log)
                valid, reason = validate_pdbqt(pdbqt)
                if not valid:
                    return row, pdbqt, log, reason
                return row, pdbqt, log, None
            except Exception as exc:
                return row, pdbqt, log, str(exc)

        if rows:
            for index, (_, future) in enumerate(bounded_futures(prepare, rows, workers,
                    cancelled=self.cancelled, name="meeko"), 1):
                row, pdbqt, log, error = future.result()
                self._emit("item_started", "meeko", row["state_id"], index, len(rows))
                try:
                    if error:
                        raise RuntimeError(error)
                    artifact_id = self.repository.add_artifact(pdbqt, "ligand_pdbqt", "meeko")
                    log_artifact = self.repository.add_artifact(log, "meeko_log", "meeko")
                    with self.repository.connection() as conn:
                        conn.execute("UPDATE conformers SET status=? WHERE conformer_id=?", ("pdbqt_ready", row["conformer_id"]))
                        conn.execute("UPDATE molecular_states SET status=?, reason=? WHERE state_id=?", ("pdbqt_ready", artifact_id, row["state_id"]))
                    prepared += 1
                    self._emit("item_succeeded", "meeko", row["state_id"], index, len(rows), succeeded=prepared)
                except Exception as exc:
                    if isinstance(exc, InterruptedError):
                        raise
                    failed += 1
                    with self.repository.connection() as conn:
                        conn.execute("UPDATE conformers SET status=? WHERE conformer_id=?", ("pdbqt_failed", row["conformer_id"]))
                        conn.execute("UPDATE molecular_states SET status=?, reason=? WHERE state_id=?", ("pdbqt_failed", str(exc)[:500], row["state_id"]))
                    self._emit("item_failed", "meeko", row["state_id"], index, len(rows), failed=failed, message=str(exc))
        self._emit("stage_completed", "meeko", total=len(rows), succeeded=prepared, skipped=skipped, failed=failed)
        return prepared, failed

    def validate_docking_configuration(self, profile_ids: Collection[str] | None = None):
        profiles = [profile for profile in self.repository.get_receptor_profiles()
                    if profile.get("enabled", False) and (profile_ids is None or str(profile.get("id")) in profile_ids)]
        if not profiles:
            raise RuntimeError("No enabled receptor profiles are configured. Open Receptors and enable an existing profile. "
                               "If the list is unexpectedly empty, check project.yml.bak and the saved receptor files before re-importing.")
        return profiles

    @exclusive_project
    def run_vina(self, profile_ids: Collection[str] | None = None) -> tuple[int, int]:
        """Dock every Meeko-prepared state against each selected receptor profile."""
        self._check_cancelled()
        profiles = self.validate_docking_configuration(profile_ids)
        # Vina is an application-level tool, shared by projects. A project-local
        # copy remains supported for portable runs.
        docking_settings = self.repository.get_settings().get("docking", {})
        backend_name = docking_settings.get("backend", "vina")
        if backend_name not in ("vina", "vinagpu"):
            raise ValueError(f"Unknown docking backend: {backend_name}")
        if backend_name == "vinagpu":
            executable = find_vinagpu_executable(self.repository.root)
            tool_hash = vinagpu_tool_hash(executable)
            backend = VinaGpuDockingBackend()
        else:
            executable = find_vina_executable(self.repository.root)
            tool_hash = file_sha256(executable)
            backend = VinaDockingBackend()
        executable_hash = file_sha256(executable)
        with self.repository.connection() as conn:
            states = conn.execute("""SELECT s.state_id, s.parent_id, a.relative_path
                FROM molecular_states s JOIN conformers c ON c.state_id=s.state_id
                JOIN artifacts a ON a.artifact_id=c.sdf_artifact_id
                WHERE s.active=1 AND c.status='pdbqt_ready' AND a.active=1 ORDER BY s.parent_id, s.state_id""").fetchall()
        campaign_settings = {
            "receptor_profiles": [{key: value for key, value in profile.items()
                                   if key not in {"reference_ligand"}} for profile in profiles],
            "receptor_profile_ids": [str(profile["id"]) for profile in profiles],
            "prepared_state_count": len(states),
            "requested_dockings": len(states) * len(profiles),
            "vina_executable": str(executable),
            "vina_executable_sha256": executable_hash,
            "docking_backend": backend_name,
            "docking_tool_hash": tool_hash,
        }
        workflow_run_id = self.repository.start_workflow_run("docking", campaign_settings)
        succeeded = failed = reused = 0
        total = len(states) * len(profiles)
        self.repository.update_workflow_stage(workflow_run_id, "vina", StageStatus.RUNNING.value,
                                              summary={"total": total})
        self.repository.record_provenance_event(
            workflow_run_id=workflow_run_id, event_type="docking_campaign_started",
            stage_name="vina", data={"total": total, "receptors": len(profiles), "states": len(states)},
        )
        self._emit("stage_started", "vina", total=total)
        index = 0
        for profile in profiles:
            self._check_cancelled()
            profile_id = str(profile["id"])
            profile_name = str(profile.get("name", profile_id))
            receptor = self.repository.root / str(profile.get("receptor", ""))
            run_settings = {key: value for key, value in profile.items()
                            if key not in {"id", "name", "enabled", "archived", "receptor"} and value is not None}
            run_settings["cpu"] = run_settings.pop("cpu_count", logical_cpu_count())
            if backend_name == "vinagpu":
                run_settings.pop("cpu", None)
                run_settings.pop("exhaustiveness", None)
                run_settings.update(backend="vinagpu", thread=docking_settings.get("thread", 8000),
                                    search_depth=docking_settings.get("search_depth", 32))
            if not receptor.is_file():
                reason = f"Prepared receptor not found: {receptor}"
                missing_fp = fingerprint(settings=run_settings, inputs={"receptor": "missing"}, tool_version=tool_hash)
                for state in states:
                    index += 1
                    failed += 1
                    run_id = str(uuid.uuid4())
                    with self.repository.connection() as conn:
                        conn.execute("UPDATE docking_runs SET is_current=0 WHERE state_id=? AND receptor_profile_id=?", (state["state_id"], profile_id))
                        conn.execute("""INSERT INTO docking_runs
                            (run_id,state_id,receptor_hash,settings_fingerprint,status,command_json,started_at,ended_at,
                             is_current,receptor_profile_id,receptor_profile_name,reason,workflow_run_id)
                            VALUES (?, ?, '', ?, ?, ?, ?, ?, 1, ?, ?, ?, ?)""",
                            (run_id, state["state_id"], missing_fp, StageStatus.FAILED,
                             json.dumps({"settings": run_settings, "argv": []}, sort_keys=True),
                             utc_now(), utc_now(), profile_id, profile_name, reason, workflow_run_id))
                    self._emit("item_failed", "vina", state["state_id"], index, total, failed=failed,
                               receptor_profile_id=profile_id, receptor_profile_name=profile_name,
                               message=reason)
                continue
            receptor_hash = file_sha256(receptor)
            parameters = _docking_parameters(run_settings, backend_name)
            for state, prepared_artifact, candidates in self._docking_plan_rows(states, profile_id):
                self._check_cancelled()
                index += 1
                self._emit("item_started", "vina", state["state_id"], index, total,
                           receptor_profile_id=profile_id, receptor_profile_name=profile_name,
                           message=f"{profile_name}: {state['state_id']}")
                ligand = self.repository.root / "artifacts" / "pdbqt" / state["parent_id"] / state["state_id"] / "ligand.pdbqt"
                run_id = str(uuid.uuid4())
                run_dir = self.repository.root / "artifacts" / "docking" / profile_id / state["parent_id"] / state["state_id"] / run_id
                output = run_dir / "vina_output.pdbqt"
                log = run_dir / "vina.log.txt"
                input_hashes = {"receptor": receptor_hash,
                                "ligand": file_sha256(ligand) if ligand.is_file() else "missing"}
                settings_fp = fingerprint(settings=parameters, inputs=input_hashes, tool_version=tool_hash)
                reusable = (self._reusable_docking(state["state_id"], profile_id, settings_fp,
                                                   parameters, receptor_hash, tool_hash, backend_name,
                                                   ligand, input_hashes["ligand"], candidates=candidates,
                                                   prepared=prepared_artifact)
                            if ligand.is_file() else None)
                if reusable:
                    if not reusable["is_current"]:
                        with self.repository.connection() as conn:
                            conn.execute("UPDATE docking_runs SET is_current=0 WHERE state_id=? AND receptor_profile_id=?",
                                         (state["state_id"], profile_id))
                            conn.execute("UPDATE docking_runs SET is_current=1 WHERE run_id=?", (reusable["run_id"],))
                    succeeded += 1
                    reused += 1
                    self._emit("item_skipped", "vina", state["state_id"], index, total, skipped=reused,
                               receptor_profile_id=profile_id, receptor_profile_name=profile_name,
                               message=f"{profile_name}: matching docking run reused")
                    continue
                try:
                    with self.repository.connection() as conn:
                        conn.execute("""INSERT INTO docking_runs
                            (run_id,state_id,receptor_hash,settings_fingerprint,status,command_json,started_at,
                             receptor_profile_id,receptor_profile_name,workflow_run_id)
                            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                            (run_id, state["state_id"], receptor_hash, settings_fp, StageStatus.RUNNING,
                             json.dumps({"settings": run_settings, "inputs": input_hashes, "argv": []}, sort_keys=True),
                             utc_now(), profile_id, profile_name, workflow_run_id))
                    result = backend.run(executable=executable, receptor=receptor, ligand=ligand, output=output, log=log, settings=run_settings, cancelled=self.cancelled)
                    raw_artifact = self.repository.add_artifact(output, "vina_output_pdbqt", "vina") if output.exists() else None
                    log_artifact = self.repository.add_artifact(log, "vina_log_txt", "vina")
                    with self.repository.connection() as conn:
                        conn.execute("UPDATE docking_runs SET is_current=0 WHERE state_id=? AND receptor_profile_id=? AND run_id<>?", (state["state_id"], profile_id, run_id))
                        conn.execute("""UPDATE docking_runs SET status=?,log_artifact_id=?,raw_output_artifact_id=?,
                            return_code=?,ended_at=?,is_current=1,reason=?,command_json=? WHERE run_id=?""",
                            (StageStatus.COMPLETED if result.poses else StageStatus.FAILED, log_artifact, raw_artifact,
                             result.return_code, utc_now(), None if result.poses else "No valid Vina poses",
                             json.dumps({"settings": run_settings, "inputs": input_hashes,
                                         "argv": list(getattr(result, "command", ()))}, sort_keys=True), run_id))
                        best_mode = min(result.poses, key=lambda pose: pose[1])[0] if result.poses else None
                        for mode, affinity, rmsd_lb, rmsd_ub in result.poses:
                            conn.execute("INSERT INTO docking_poses(pose_id,run_id,mode_index,affinity,rmsd_lb,rmsd_ub,is_best_for_state) VALUES (?, ?, ?, ?, ?, ?, ?)", (str(uuid.uuid4()), run_id, mode, affinity, rmsd_lb, rmsd_ub, int(mode == best_mode)))
                    succeeded += 1 if result.poses else 0
                    failed += 0 if result.poses else 1
                    self._emit("item_succeeded" if result.poses else "item_failed", "vina", state["state_id"], index, total,
                               succeeded=succeeded, failed=failed, receptor_profile_id=profile_id, receptor_profile_name=profile_name,
                               message=f"{profile_name}: complete" if result.poses else f"{profile_name}: no valid Vina poses")
                except Exception as exc:
                    if isinstance(exc, InterruptedError):
                        with self.repository.connection() as conn:
                            conn.execute("UPDATE docking_runs SET status=?, ended_at=?, reason=?, is_current=0 WHERE run_id=?",
                                         (StageStatus.CANCELLED, utc_now(), "Pipeline stopped safely", run_id))
                        self.repository.update_workflow_stage(workflow_run_id, "vina", StageStatus.CANCELLED.value,
                                                              summary={"total": total, "succeeded": succeeded, "failed": failed, "reused": reused},
                                                              error="Pipeline stopped safely")
                        self.repository.finish_workflow_run(workflow_run_id, StageStatus.CANCELLED.value,
                                                            error="Pipeline stopped safely")
                        raise
                    failed += 1
                    with self.repository.connection() as conn:
                        conn.execute("UPDATE docking_runs SET is_current=0 WHERE state_id=? AND receptor_profile_id=? AND run_id<>?", (state["state_id"], profile_id, run_id))
                        conn.execute("UPDATE docking_runs SET status=?, ended_at=?, reason=?, is_current=1 WHERE run_id=?", (StageStatus.FAILED, utc_now(), str(exc)[:500], run_id))
                    self._emit("item_failed", "vina", state["state_id"], index, total, failed=failed,
                               receptor_profile_id=profile_id, receptor_profile_name=profile_name,
                               message=f"{profile_name}: {exc}")
        self._emit("stage_completed", "vina", total=total, succeeded=succeeded, skipped=reused, failed=failed)
        summary = {"total": total, "succeeded": succeeded, "failed": failed, "reused": reused,
                   "newly_completed": max(0, succeeded - reused), "receptor_count": len(profiles),
                   "state_count": len(states)}
        final_status = StageStatus.COMPLETED.value if failed == 0 else StageStatus.FAILED.value
        self.repository.update_workflow_stage(workflow_run_id, "vina", final_status, summary=summary,
                                              error=None if failed == 0 else f"{failed} docking item(s) failed")
        self.repository.record_provenance_event(
            workflow_run_id=workflow_run_id, event_type="docking_campaign_finished",
            stage_name="vina", level="INFO" if failed == 0 else "WARNING", data=summary,
        )
        self.repository.finish_workflow_run(workflow_run_id, final_status,
                                            error=None if failed == 0 else f"{failed} docking item(s) failed")
        self.repository.record_tool(backend_name, executable.name, executable)
        return succeeded, failed

    @exclusive_project
    def run_postdock(self) -> tuple[int, int]:
        self._check_cancelled()
        settings = self.repository.get_settings().get("postdock", {})
        mode = str(settings.get("mode", "split_and_sdf"))
        limit = max(1, int(settings.get("poses_per_compound", 3)))
        selected = set(settings.get("selected_parents", []))
        profiles = [profile for profile in self.repository.get_receptor_profiles() if profile.get("enabled")]
        profiles_by_id = {str(profile["id"]): profile for profile in profiles}
        profile_ids = list(profiles_by_id)
        split_executable = find_vina_split_executable(self.repository.root) if mode in ("split_only", "split_and_sdf") else None
        with self.repository.connection() as conn:
            placeholders = ",".join("?" for _ in profile_ids)
            runs = conn.execute(f"""SELECT d.run_id, d.receptor_profile_id, d.receptor_profile_name,
                    s.parent_id, s.state_id, a.relative_path raw_path
                FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id
                JOIN artifacts a ON a.artifact_id=d.raw_output_artifact_id
                WHERE d.status='completed' AND d.is_current=1 AND s.active=1
                AND d.receptor_profile_id IN ({placeholders}) ORDER BY d.receptor_profile_id, s.parent_id, s.state_id""", profile_ids).fetchall() if profile_ids else []
        if selected:
            runs = [run for run in runs if run["parent_id"] in selected]
        self._emit("stage_started", "postdock", total=len(runs))
        service = PostDockService()
        for profile_id in sorted({str(run["receptor_profile_id"]) for run in runs}):
            self._check_cancelled()
            profile = profiles_by_id[profile_id]
            receptor_pdbqt = self.repository.root / str(profile["receptor"])
            prepared_pdb = receptor_pdbqt.parent / "receptor_prepared.pdb"
            profile_root = self.repository.root / "For_PostDocking" / profile_id
            receptor_export = profile_root / receptor_export_filename(str(profile.get("name", profile_id)))
            source_kind = service.export_receptor_pdb(receptor_pdbqt, receptor_export, prepared_pdb)
            self.repository.add_artifact(receptor_export, "postdock_receptor_pdb", "postdock")
            self.repository.record_provenance_event(
                event_type="postdock_receptor_exported", stage_name="postdock",
                entity_type="receptor_profile", entity_id=profile_id, receptor_profile_id=profile_id,
                data={"source_kind": source_kind, "source": str(prepared_pdb if source_kind == "meeko_prepared_pdb" else receptor_pdbqt),
                      "output": receptor_export.relative_to(self.repository.root).as_posix()},
            )
            with self.repository.connection() as conn:
                validation_runs = conn.execute("""SELECT run_id FROM redocking_runs
                    WHERE receptor_profile_id=? AND status='ARTIFACTS_READY'
                    ORDER BY COALESCE(finished_at,started_at) DESC,run_id DESC""", (profile_id,)).fetchall()
            standard_source = next((self.repository.root / "inputs" / "receptors" / profile_id / "redocking"
                                    / str(row["run_id"]) / "dockrmsd" / "reference_ligand_heavy.mol2"
                                    for row in validation_runs
                                    if (self.repository.root / "inputs" / "receptors" / profile_id / "redocking"
                                        / str(row["run_id"]) / "dockrmsd" / "reference_ligand_heavy.mol2").is_file()), None)
            if standard_source:
                standard_export = profile_root / "reference_ligand_heavy.mol2"
                service.export_validation_standard(standard_source, standard_export)
                self.repository.add_artifact(standard_export, "postdock_validation_standard_mol2", "postdock")
                self.repository.record_provenance_event(
                    event_type="postdock_validation_standard_exported", stage_name="postdock",
                    entity_type="receptor_profile", entity_id=profile_id, receptor_profile_id=profile_id,
                    data={"source": standard_source.relative_to(self.repository.root).as_posix(),
                          "output": standard_export.relative_to(self.repository.root).as_posix()},
                )
            else:
                self.repository.record_provenance_event(
                    event_type="postdock_validation_standard_unavailable", stage_name="postdock", level="WARNING",
                    entity_type="receptor_profile", entity_id=profile_id, receptor_profile_id=profile_id,
                    message="No completed receptor-validation reference MOL2 was available for export",
                )
        success = failed = 0
        for run in runs:
            self._check_cancelled()
            profile_root = self.repository.root / "For_PostDocking" / run["receptor_profile_id"]
            expected_folders = ("PDBQTs",) if mode == "split_only" else ("SDF",) if mode == "sdf_only" else ("SDF", "PDBQTs")
            export_root = lambda folder: profile_root / folder / run["parent_id"] / run["state_id"] / run["run_id"]
            already_exported = all(
                export_root(folder).is_dir()
                and any(export_root(folder).glob("*.sdf" if folder == "SDF" else "*.pdbqt"))
                for folder in expected_folders
            )
            if already_exported:
                self._emit("item_skipped", "postdock", run["state_id"], success + failed, len(runs), skipped=1, message="Compound already exported")
                continue
            raw = self.repository.root / run["raw_path"]
            base_dir = profile_root
            pdbqt_dir = base_dir / "PDBQTs" / run["parent_id"] / run["state_id"] / run["run_id"]
            sdf_dir = base_dir / "SDF" / run["parent_id"] / run["state_id"] / run["run_id"]
            try:
                pose_files = []
                if mode in ("split_only", "split_and_sdf"):
                    pose_files = service.split(split_executable, raw, pdbqt_dir, f"{run['parent_id']}_{run['state_id']}_pose_")[:limit]
                    if not pose_files:
                        raise RuntimeError("vina_split produced no pose files")
                if mode in ("split_and_sdf", "sdf_only"):
                    sources = pose_files if pose_files else [raw]
                    for source in sources:
                        suffix = source.stem if len(sources) > 1 else "all"
                        service.export_sdf(source, sdf_dir / f"{suffix}.sdf", sdf_dir / f"{suffix}.export.log.txt")
                success += 1
                self._emit("item_succeeded", "postdock", run["state_id"], success, len(runs), succeeded=success)
            except Exception as exc:
                failed += 1
                error_dir = sdf_dir if mode in ("split_and_sdf", "sdf_only") else pdbqt_dir
                error_dir.mkdir(parents=True, exist_ok=True)
                (error_dir / "postdock.error.txt").write_text(str(exc), encoding="utf-8")
                self._emit("item_failed", "postdock", run["state_id"], success + failed, len(runs), failed=failed, message=str(exc))
        self._emit("stage_completed", "postdock", total=len(runs), succeeded=success, failed=failed)
        return success, failed
