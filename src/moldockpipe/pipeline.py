from __future__ import annotations

import json
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Collection

from .fingerprints import file_sha256, fingerprint
from .models import ProgressCallback, ProgressEvent, ScreeningPolicy, StageStatus
from .project import ProjectRepository, utc_now
from .services.molscrub import MolScrubService
from .services.meeko import MeekoService
from .services.screening import screen_smiles
from .services.validation import validate_pdbqt
from .services.vina import VinaDockingBackend
from .services.postdock import PostDockService


class PipelineRunner:
    """Application service for small, independently resumable stage operations."""

    def __init__(self, repository: ProjectRepository, progress: ProgressCallback | None = None):
        self.repository = repository
        self.progress = progress

    def _emit(self, event: str, stage: str, item_id: str | None = None, index: int = 0, total: int = 0, **kwargs: object) -> None:
        if self.progress:
            self.progress(ProgressEvent(event, stage, item_id, index, total, **kwargs))

    def run_screening(self) -> tuple[int, int]:
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
            parents = conn.execute("SELECT parent_id, source_smiles FROM parent_ligands WHERE active=1 ORDER BY parent_id").fetchall()
        self._emit("stage_started", "screening", total=len(parents))
        skipped = 0
        for index, parent in enumerate(parents, 1):
            self._emit("item_started", "screening", parent["parent_id"], index, len(parents))
            result = screen_smiles(parent["source_smiles"], enabled, policy)
            run_fingerprint = fingerprint(settings={"enabled": enabled, "policy": policy.value}, inputs={"smiles": parent["source_smiles"]}, tool_version=rdkit_version)
            with self.repository.connection() as conn:
                previous = conn.execute("SELECT fingerprint,status FROM screening_results WHERE parent_id=? AND active=1", (parent["parent_id"],)).fetchone()
            if previous and previous["fingerprint"] == run_fingerprint and previous["status"] == StageStatus.COMPLETED:
                skipped += 1
                completed += 1
                self._emit("item_skipped", "screening", parent["parent_id"], index, len(parents), skipped=skipped, message="Matching screening result reused")
                continue
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

    def run_molscrub(self) -> tuple[int, int]:
        project = self.repository.get_settings()
        settings = project.get("molscrub", {})
        screening = project.get("screening", {})
        policy = screening.get("policy", ScreeningPolicy.ANNOTATE_ONLY.value)
        with self.repository.connection() as conn:
            rows = conn.execute("""SELECT p.parent_id, p.source_smiles, s.decision FROM parent_ligands p
                LEFT JOIN screening_results s ON s.parent_id=p.parent_id AND s.active=1
                WHERE p.active=1 ORDER BY p.parent_id""").fetchall()
        allowed = [row for row in rows if row["decision"] not in ("fail", None) or policy == ScreeningPolicy.ANNOTATE_ONLY.value]
        service = MolScrubService()
        created = failed = skipped = 0
        self._emit("stage_started", "molscrub", total=len(allowed))
        for index, row in enumerate(allowed, 1):
            self._emit("item_started", "molscrub", row["parent_id"], index, len(allowed))
            with self.repository.connection() as conn:
                existing = conn.execute("""SELECT COUNT(*) FROM molecular_states s JOIN conformers c ON c.state_id=s.state_id
                    JOIN artifacts a ON a.artifact_id=c.sdf_artifact_id WHERE s.parent_id=? AND s.active=1 AND c.status IN ('ready','pdbqt_ready') AND a.active=1""", (row["parent_id"],)).fetchone()[0]
            if existing:
                created += existing
                skipped += 1
                self._emit("item_skipped", "molscrub", row["parent_id"], index, len(allowed), skipped=skipped, message=f"{existing} existing states reused")
                continue
            try:
                states, truncated = service.generate(row["source_smiles"], ph=float(settings.get("ph", 7.4)), enumerate_states=bool(settings.get("enumerate_states", True)), max_states=int(settings.get("max_states", 32)))
                self._persist_states(row["parent_id"], states, truncated, settings)
                created += len(states)
                self._emit("item_succeeded", "molscrub", row["parent_id"], index, len(allowed), succeeded=1)
            except Exception as exc:
                failed += 1
                with self.repository.connection() as conn:
                    conn.execute("UPDATE parent_ligands SET parse_reason=? WHERE parent_id=?", (f"MolScrub failed: {exc}", row["parent_id"]))
                self._emit("item_failed", "molscrub", row["parent_id"], index, len(allowed), failed=1, message=str(exc))
        self._emit("stage_completed", "molscrub", total=len(allowed), succeeded=created, skipped=skipped, failed=failed)
        return created, failed

    def _persist_states(self, parent_id: str, states: list[object], truncated: bool, settings: dict[str, object]) -> None:
        from rdkit import Chem
        try:
            import molscrub
            version = getattr(molscrub, "__version__", "unknown")
        except ImportError:
            version = "unknown"
        for state in states:
            state_id = f"{parent_id}-{state.structure_hash[:16]}"
            state_fingerprint = fingerprint(settings=settings, inputs={"state": state.structure_hash}, tool_version=version)
            state_dir = self.repository.root / "artifacts" / "sdf" / parent_id / state_id
            state_dir.mkdir(parents=True, exist_ok=True)
            sdf_path = state_dir / "input_conformer_0001.sdf"
            writer = Chem.SDWriter(str(sdf_path))
            writer.write(state.molecule)
            writer.close()
            with self.repository.connection() as conn:
                conn.execute("""INSERT OR REPLACE INTO molecular_states
                    (state_id,parent_id,state_smiles,state_isomeric_smiles,state_inchikey,formal_charge,tautomer_index,protomer_index,state_structure_hash,generation_fingerprint,status,reason,enumeration_truncated,created_at)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", (state_id, parent_id, state.state_smiles, state.state_isomeric_smiles, state.state_inchikey, state.formal_charge, state.tautomer_index, state.protomer_index, state.structure_hash, state_fingerprint, "sdf_ready", None, int(truncated), utc_now()))
            artifact_id = self.repository.add_artifact(sdf_path, "state_sdf", "molscrub")
            with self.repository.connection() as conn:
                conn.execute("INSERT OR REPLACE INTO conformers VALUES (?, ?, ?, ?, ?, ?)", (f"{state_id}-conf-1", state_id, 1, artifact_id, "ready", utc_now()))

    def run_meeko(self) -> tuple[int, int]:
        """Prepare every retained SDF state into a state-specific ligand PDBQT."""
        prepared = failed = 0
        settings = self.repository.get_settings().get("meeko", {})
        workers = max(1, min(32, int(settings.get("workers", 4))))
        with self.repository.connection() as conn:
            all_states = conn.execute("SELECT COUNT(*) FROM molecular_states WHERE active=1").fetchone()[0]
            rows = conn.execute("""SELECT s.state_id, s.parent_id, c.conformer_id, a.relative_path
                FROM molecular_states s JOIN conformers c ON c.state_id=s.state_id
                JOIN artifacts a ON a.artifact_id=c.sdf_artifact_id WHERE s.active=1 AND c.status='ready' AND a.active=1
                ORDER BY s.parent_id, s.state_id""").fetchall()
        self._emit("stage_started", "meeko", total=len(rows))
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

        with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="meeko") as pool:
            futures = [pool.submit(prepare, row) for row in rows]
            for index, future in enumerate(as_completed(futures), 1):
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
                    failed += 1
                    with self.repository.connection() as conn:
                        conn.execute("UPDATE conformers SET status=? WHERE conformer_id=?", ("pdbqt_failed", row["conformer_id"]))
                        conn.execute("UPDATE molecular_states SET status=?, reason=? WHERE state_id=?", ("pdbqt_failed", str(exc)[:500], row["state_id"]))
                    self._emit("item_failed", "meeko", row["state_id"], index, len(rows), failed=failed, message=str(exc))
        self._emit("stage_completed", "meeko", total=len(rows), succeeded=prepared, skipped=skipped, failed=failed)
        return prepared, failed

    def run_vina(self, profile_ids: Collection[str] | None = None) -> tuple[int, int]:
        """Dock every Meeko-prepared state against each selected receptor profile."""
        profiles = [profile for profile in self.repository.get_receptor_profiles()
                    if profile.get("enabled", False) and (profile_ids is None or str(profile.get("id")) in profile_ids)]
        if not profiles:
            raise RuntimeError("No enabled receptor profiles are configured")
        # Vina is an application-level tool, shared by projects. A project-local
        # copy remains supported for portable runs.
        application_root = Path(__file__).resolve().parents[2]
        candidates = [
            application_root / "tools" / "vina" / "vina.exe",
            application_root / "tools" / "vina" / "vina_1.2.7_win.exe",
            self.repository.root / "tools" / "vina" / "vina.exe",
            self.repository.root / "tools" / "vina" / "vina_1.2.7_win.exe",
        ]
        executable = next((candidate for candidate in candidates if candidate.is_file()), candidates[0])
        if not executable.is_file():
            raise FileNotFoundError(f"Vina executable not found. Place it at {application_root / 'tools' / 'vina'}")
        executable_hash = file_sha256(executable)
        with self.repository.connection() as conn:
            states = conn.execute("""SELECT s.state_id, s.parent_id, a.relative_path
                FROM molecular_states s JOIN conformers c ON c.state_id=s.state_id
                JOIN artifacts a ON a.artifact_id=c.sdf_artifact_id
                WHERE s.active=1 AND c.status='pdbqt_ready' AND a.active=1 ORDER BY s.parent_id, s.state_id""").fetchall()
        backend = VinaDockingBackend()
        succeeded = failed = 0
        total = len(states) * len(profiles)
        self._emit("stage_started", "vina", total=total)
        index = 0
        for profile in profiles:
            profile_id = str(profile["id"])
            profile_name = str(profile.get("name", profile_id))
            receptor = self.repository.root / str(profile.get("receptor", ""))
            run_settings = {key: value for key, value in profile.items()
                            if key not in {"id", "name", "enabled", "archived", "receptor"} and value is not None}
            run_settings["cpu"] = run_settings.pop("cpu_count", 1)
            if not receptor.is_file():
                reason = f"Prepared receptor not found: {receptor}"
                missing_fp = fingerprint(settings=run_settings, inputs={"receptor": "missing"}, tool_version=executable_hash)
                for state in states:
                    index += 1
                    failed += 1
                    run_id = str(uuid.uuid4())
                    with self.repository.connection() as conn:
                        conn.execute("UPDATE docking_runs SET is_current=0 WHERE state_id=? AND receptor_profile_id=?", (state["state_id"], profile_id))
                        conn.execute("""INSERT INTO docking_runs
                            (run_id,state_id,receptor_hash,settings_fingerprint,status,command_json,started_at,ended_at,
                             is_current,receptor_profile_id,receptor_profile_name,reason)
                            VALUES (?, ?, '', ?, ?, ?, ?, ?, 1, ?, ?, ?)""",
                            (run_id, state["state_id"], missing_fp, StageStatus.FAILED, json.dumps(run_settings), utc_now(), utc_now(), profile_id, profile_name, reason))
                    self._emit("item_failed", "vina", state["state_id"], index, total, failed=failed,
                               receptor_profile_id=profile_id, receptor_profile_name=profile_name,
                               message=reason)
                continue
            receptor_hash = file_sha256(receptor)
            settings_fp = fingerprint(settings=run_settings, inputs={"receptor": receptor_hash}, tool_version=executable_hash)
            for state in states:
                index += 1
                self._emit("item_started", "vina", state["state_id"], index, total,
                           receptor_profile_id=profile_id, receptor_profile_name=profile_name,
                           message=f"{profile_name}: {state['state_id']}")
                ligand = self.repository.root / "artifacts" / "pdbqt" / state["parent_id"] / state["state_id"] / "ligand.pdbqt"
                run_id = str(uuid.uuid4())
                run_dir = self.repository.root / "artifacts" / "docking" / profile_id / state["parent_id"] / state["state_id"] / run_id
                output = run_dir / "vina_output.pdbqt"
                log = run_dir / "vina.log.txt"
                with self.repository.connection() as conn:
                    reusable = conn.execute("""SELECT d.run_id, raw.relative_path AS raw_path, log.relative_path AS log_path
                        FROM docking_runs d LEFT JOIN artifacts raw ON raw.artifact_id=d.raw_output_artifact_id
                        LEFT JOIN artifacts log ON log.artifact_id=d.log_artifact_id
                        WHERE d.state_id=? AND d.receptor_profile_id=? AND d.settings_fingerprint=? AND d.status=? AND d.is_current=1""",
                        (state["state_id"], profile_id, settings_fp, StageStatus.COMPLETED)).fetchone()
                if reusable and reusable["raw_path"] and reusable["log_path"] and (self.repository.root / reusable["raw_path"]).is_file() and (self.repository.root / reusable["log_path"]).is_file():
                    succeeded += 1
                    self._emit("item_skipped", "vina", state["state_id"], index, total, skipped=1,
                               receptor_profile_id=profile_id, receptor_profile_name=profile_name,
                               message=f"{profile_name}: matching docking run reused")
                    continue
                try:
                    with self.repository.connection() as conn:
                        conn.execute("""INSERT INTO docking_runs
                            (run_id,state_id,receptor_hash,settings_fingerprint,status,command_json,started_at,receptor_profile_id,receptor_profile_name)
                            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                            (run_id, state["state_id"], receptor_hash, settings_fp, StageStatus.RUNNING,
                             json.dumps(run_settings), utc_now(), profile_id, profile_name))
                    result = backend.run(executable=executable, receptor=receptor, ligand=ligand, output=output, log=log, settings=run_settings)
                    raw_artifact = self.repository.add_artifact(output, "vina_output_pdbqt", "vina") if output.exists() else None
                    log_artifact = self.repository.add_artifact(log, "vina_log_txt", "vina")
                    with self.repository.connection() as conn:
                        conn.execute("UPDATE docking_runs SET is_current=0 WHERE state_id=? AND receptor_profile_id=? AND run_id<>?", (state["state_id"], profile_id, run_id))
                        conn.execute("UPDATE docking_runs SET status=?, log_artifact_id=?, raw_output_artifact_id=?, return_code=?, ended_at=?, is_current=1, reason=? WHERE run_id=?", (StageStatus.COMPLETED if result.poses else StageStatus.FAILED, log_artifact, raw_artifact, result.return_code, utc_now(), None if result.poses else "No valid Vina poses", run_id))
                        best_mode = min(result.poses, key=lambda pose: pose[1])[0] if result.poses else None
                        for mode, affinity, rmsd_lb, rmsd_ub in result.poses:
                            conn.execute("INSERT INTO docking_poses(pose_id,run_id,mode_index,affinity,rmsd_lb,rmsd_ub,is_best_for_state) VALUES (?, ?, ?, ?, ?, ?, ?)", (str(uuid.uuid4()), run_id, mode, affinity, rmsd_lb, rmsd_ub, int(mode == best_mode)))
                    succeeded += 1 if result.poses else 0
                    failed += 0 if result.poses else 1
                    self._emit("item_succeeded" if result.poses else "item_failed", "vina", state["state_id"], index, total,
                               succeeded=succeeded, failed=failed, receptor_profile_id=profile_id, receptor_profile_name=profile_name,
                               message=f"{profile_name}: complete" if result.poses else f"{profile_name}: no valid Vina poses")
                except Exception as exc:
                    failed += 1
                    with self.repository.connection() as conn:
                        conn.execute("UPDATE docking_runs SET is_current=0 WHERE state_id=? AND receptor_profile_id=? AND run_id<>?", (state["state_id"], profile_id, run_id))
                        conn.execute("UPDATE docking_runs SET status=?, ended_at=?, reason=?, is_current=1 WHERE run_id=?", (StageStatus.FAILED, utc_now(), str(exc)[:500], run_id))
                    self._emit("item_failed", "vina", state["state_id"], index, total, failed=failed,
                               receptor_profile_id=profile_id, receptor_profile_name=profile_name,
                               message=f"{profile_name}: {exc}")
        self._emit("stage_completed", "vina", total=total, succeeded=succeeded, failed=failed)
        self.repository.record_tool("vina", executable.name, executable)
        return succeeded, failed

    def run_postdock(self) -> tuple[int, int]:
        settings = self.repository.get_settings().get("postdock", {})
        mode = str(settings.get("mode", "split_and_sdf"))
        limit = max(1, int(settings.get("poses_per_compound", 3)))
        selected = set(settings.get("selected_parents", []))
        profile_ids = [str(profile["id"]) for profile in self.repository.get_receptor_profiles() if profile.get("enabled")]
        application_root = Path(__file__).resolve().parents[2]
        split_executable = next((p for p in (application_root / "tools" / "vina" / "vina_split.exe", application_root / "tools" / "vina_split.exe") if p.is_file()), None)
        if mode in ("split_only", "split_and_sdf") and split_executable is None:
            raise FileNotFoundError("vina_split.exe not found in tools/vina or tools")
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
        success = failed = 0
        for run in runs:
            profile_root = self.repository.root / "For_PostDocking" / run["receptor_profile_id"]
            expected_folders = ("PDBQTs",) if mode == "split_only" else ("SDF",) if mode == "sdf_only" else ("SDF", "PDBQTs")
            already_exported = all((profile_root / folder / run["parent_id"]).exists() for folder in expected_folders)
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
