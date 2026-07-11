from __future__ import annotations

import json
import uuid
from pathlib import Path

from .fingerprints import fingerprint
from .models import ScreeningPolicy, StageStatus
from .project import ProjectRepository, utc_now
from .services.molscrub import MolScrubService
from .services.meeko import MeekoService
from .services.screening import screen_smiles
from .services.validation import validate_pdbqt
from .services.vina import VinaDockingBackend


class PipelineRunner:
    """Application service for small, independently resumable stage operations."""

    def __init__(self, repository: ProjectRepository):
        self.repository = repository

    def run_screening(self) -> tuple[int, int]:
        settings = self.repository.get_settings().get("screening", {})
        policy = ScreeningPolicy(settings.get("policy", ScreeningPolicy.ANNOTATE_ONLY.value))
        enabled = {name: bool(settings.get(name, False)) for name in ("lipinski", "veber", "egan", "ghose")}
        try:
            import rdkit
            rdkit_version = rdkit.__version__
        except ImportError as exc:
            raise RuntimeError("RDKit must be installed before screening") from exc
        completed = failed = 0
        with self.repository.connection() as conn:
            parents = conn.execute("SELECT parent_id, source_smiles FROM parent_ligands ORDER BY parent_id").fetchall()
        for parent in parents:
            result = screen_smiles(parent["source_smiles"], enabled, policy)
            run_fingerprint = fingerprint(settings={"enabled": enabled, "policy": policy.value}, inputs={"smiles": parent["source_smiles"]}, tool_version=rdkit_version)
            run_id = str(uuid.uuid4())
            with self.repository.connection() as conn:
                conn.execute("INSERT INTO stage_runs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", (run_id, "screening", "parent", parent["parent_id"], run_fingerprint, StageStatus.RUNNING, utc_now(), None, None))
                conn.execute("UPDATE parent_ligands SET canonical_source_smiles=?, parent_inchikey=?, parse_status=?, parse_reason=? WHERE parent_id=?", (result.canonical_smiles, result.inchikey, StageStatus.COMPLETED if result.decision != "fail" else StageStatus.FAILED, result.reason, parent["parent_id"]))
                conn.execute("""INSERT OR REPLACE INTO screening_results VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """, (parent["parent_id"], run_fingerprint, StageStatus.COMPLETED if result.decision != "fail" else StageStatus.FAILED, json.dumps(result.descriptors), json.dumps(result.rules), result.decision, result.reason, utc_now()))
                conn.execute("UPDATE stage_runs SET status=?, ended_at=?, reason=? WHERE stage_run_id=?", (StageStatus.COMPLETED if result.decision != "fail" else StageStatus.FAILED, utc_now(), result.reason, run_id))
            if result.decision == "fail":
                failed += 1
            else:
                completed += 1
        self.repository.record_tool("rdkit", rdkit_version, None)
        return completed, failed

    def run_molscrub(self) -> tuple[int, int]:
        project = self.repository.get_settings()
        settings = project.get("molscrub", {})
        screening = project.get("screening", {})
        policy = screening.get("policy", ScreeningPolicy.ANNOTATE_ONLY.value)
        with self.repository.connection() as conn:
            rows = conn.execute("""SELECT p.parent_id, p.source_smiles, s.decision FROM parent_ligands p
                LEFT JOIN screening_results s USING(parent_id) ORDER BY p.parent_id""").fetchall()
        allowed = [row for row in rows if row["decision"] not in ("fail", None) or policy == ScreeningPolicy.ANNOTATE_ONLY.value]
        service = MolScrubService()
        created = failed = 0
        for row in allowed:
            try:
                states, truncated = service.generate(row["source_smiles"], ph=float(settings.get("ph", 7.4)), enumerate_states=bool(settings.get("enumerate_states", True)), max_states=int(settings.get("max_states", 32)))
                self._persist_states(row["parent_id"], states, truncated, settings)
                created += len(states)
            except Exception as exc:
                failed += 1
                with self.repository.connection() as conn:
                    conn.execute("UPDATE parent_ligands SET parse_reason=? WHERE parent_id=?", (f"MolScrub failed: {exc}", row["parent_id"]))
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
        service = MeekoService()
        prepared = failed = 0
        with self.repository.connection() as conn:
            rows = conn.execute("""SELECT s.state_id, s.parent_id, c.conformer_id, a.relative_path
                FROM molecular_states s JOIN conformers c ON c.state_id=s.state_id
                JOIN artifacts a ON a.artifact_id=c.sdf_artifact_id WHERE c.status='ready'
                ORDER BY s.parent_id, s.state_id""").fetchall()
        for row in rows:
            sdf = self.repository.root / row["relative_path"]
            output_dir = self.repository.root / "artifacts" / "pdbqt" / row["parent_id"] / row["state_id"]
            pdbqt = output_dir / "ligand.pdbqt"
            log = self.repository.root / "logs" / "meeko" / f"{row['state_id']}.log"
            try:
                service.prepare_ligand(sdf, pdbqt, log)
                valid, reason = validate_pdbqt(pdbqt)
                if not valid:
                    raise RuntimeError(reason)
                artifact_id = self.repository.add_artifact(pdbqt, "ligand_pdbqt", "meeko")
                log_artifact = self.repository.add_artifact(log, "meeko_log", "meeko")
                with self.repository.connection() as conn:
                    conn.execute("UPDATE conformers SET status=? WHERE conformer_id=?", ("pdbqt_ready", row["conformer_id"]))
                    conn.execute("UPDATE molecular_states SET status=?, reason=? WHERE state_id=?", ("pdbqt_ready", artifact_id, row["state_id"]))
                prepared += 1
            except Exception as exc:
                failed += 1
                with self.repository.connection() as conn:
                    conn.execute("UPDATE conformers SET status=? WHERE conformer_id=?", ("pdbqt_failed", row["conformer_id"]))
                    conn.execute("UPDATE molecular_states SET status=?, reason=? WHERE state_id=?", ("pdbqt_failed", str(exc)[:500], row["state_id"]))
        return prepared, failed

    def run_vina(self) -> tuple[int, int]:
        """Dock every Meeko-prepared state and persist all valid Vina modes."""
        settings = dict(self.repository.get_settings().get("vina", {}))
        # Vina is an application-level tool, shared by projects. A project-local
        # copy remains supported for portable runs.
        settings.pop("executable", None)  # legacy project.yml field; discovery is fixed
        application_root = Path(__file__).resolve().parents[2]
        candidates = [
            application_root / "tools" / "vina" / "vina.exe",
            application_root / "tools" / "vina" / "vina_1.2.7_win.exe",
            self.repository.root / "tools" / "vina" / "vina.exe",
            self.repository.root / "tools" / "vina" / "vina_1.2.7_win.exe",
        ]
        executable = next((candidate for candidate in candidates if candidate.is_file()), candidates[0])
        receptor = self.repository.root / str(settings.pop("receptor", "inputs/receptor_prepared.pdbqt"))
        if not executable.is_file():
            raise FileNotFoundError(f"Vina executable not found. Place it at {application_root / 'tools' / 'vina'}")
        if not receptor.is_file():
            raise FileNotFoundError(f"Prepared receptor not found: {receptor}")
        with self.repository.connection() as conn:
            states = conn.execute("""SELECT s.state_id, s.parent_id, a.relative_path
                FROM molecular_states s JOIN conformers c ON c.state_id=s.state_id
                JOIN artifacts a ON a.artifact_id=c.sdf_artifact_id
                WHERE c.status='pdbqt_ready' ORDER BY s.parent_id, s.state_id""").fetchall()
            pdbqt_artifacts = {row["relative_path"].split("/")[-3]: row for row in conn.execute("SELECT artifact_id, relative_path FROM artifacts WHERE artifact_type='ligand_pdbqt' AND active=1")}
        backend = VinaDockingBackend()
        succeeded = failed = 0
        receptor_hash = __import__("moldockpipe.fingerprints", fromlist=["file_sha256"]).file_sha256(receptor)
        for state in states:
            # The PDBQT artifact directory is deterministic and state-specific.
            ligand = self.repository.root / "artifacts" / "pdbqt" / state["parent_id"] / state["state_id"] / "ligand.pdbqt"
            run_id = str(uuid.uuid4())
            run_dir = self.repository.root / "artifacts" / "docking" / state["parent_id"] / state["state_id"] / run_id
            output = run_dir / "vina_output.pdbqt"
            log = run_dir / "vina.log.txt"
            run_settings = {key: value for key, value in settings.items() if value is not None}
            try:
                with self.repository.connection() as conn:
                    conn.execute("INSERT INTO docking_runs(run_id,state_id,receptor_hash,settings_fingerprint,status,command_json,started_at) VALUES (?, ?, ?, ?, ?, ?, ?)", (run_id, state["state_id"], receptor_hash, fingerprint(settings=run_settings, inputs={"receptor": receptor_hash}, tool_version=str(executable)), StageStatus.RUNNING, json.dumps(run_settings), utc_now()))
                result = backend.run(executable=executable, receptor=receptor, ligand=ligand, output=output, log=log, settings=run_settings)
                raw_artifact = self.repository.add_artifact(output, "vina_output_pdbqt", "vina") if output.exists() else None
                log_artifact = self.repository.add_artifact(log, "vina_log_txt", "vina")
                with self.repository.connection() as conn:
                    conn.execute("UPDATE docking_runs SET status=?, log_artifact_id=?, raw_output_artifact_id=?, return_code=?, ended_at=?, is_current=? WHERE run_id=?", (StageStatus.COMPLETED if result.poses else StageStatus.FAILED, log_artifact, raw_artifact, result.return_code, utc_now(), int(bool(result.poses)), run_id))
                    for mode, affinity, rmsd_lb, rmsd_ub in result.poses:
                        conn.execute("INSERT INTO docking_poses(pose_id,run_id,mode_index,affinity,rmsd_lb,rmsd_ub,is_best_for_state) VALUES (?, ?, ?, ?, ?, ?, ?)", (str(uuid.uuid4()), run_id, mode, affinity, rmsd_lb, rmsd_ub, int(mode == min(result.poses, key=lambda pose: pose[1])[0])))
                succeeded += 1 if result.poses else 0
                failed += 0 if result.poses else 1
            except Exception as exc:
                failed += 1
                with self.repository.connection() as conn:
                    conn.execute("UPDATE docking_runs SET status=?, ended_at=? WHERE run_id=?", (StageStatus.FAILED, utc_now(), run_id))
        return succeeded, failed
