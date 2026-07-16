from __future__ import annotations

import csv
import json
import shutil
import threading
import uuid
from importlib.metadata import version
from pathlib import Path
from typing import Any, Callable

from ..fingerprints import file_sha256, fingerprint
from ..project import ProjectRepository, utc_now
from ..receptors.ligand_chemistry import (
    chemistry_summary, hydrogenated_copy, read_sdf, validate_same_heavy_graph, write_mol2, write_sdf,
)
from ..services.meeko import MeekoService
from ..services.postdock import PostDockService
from ..services.validation import parse_vina_poses, validate_pdbqt
from ..services.vina import VinaDockingBackend, find_vina_executable
from .models import RedockingSettings, RedockingStatus

STAGES = ("reference_validation", "ligand_preparation", "vina_redocking", "pose_export", "mol2_generation", "final_validation")


def _tool_version(name: str) -> str:
    try: return version(name)
    except Exception: return "unknown"


def _profile_path(repository: ProjectRepository, text: object) -> Path:
    path = Path(str(text or "")); return path if path.is_absolute() else repository.root / path


def validate_redocking_prerequisites(repository: ProjectRepository, profile: dict[str, Any]) -> list[str]:
    missing: list[str] = []
    receptor = _profile_path(repository, profile.get("receptor"))
    valid, _ = validate_pdbqt(receptor)
    if not valid: missing.append("Valid prepared receptor PDBQT")
    reference = profile.get("reference_ligand")
    if not isinstance(reference, dict) or not reference.get("identity"): missing.append("Reference-ligand identity")
    sdf = _profile_path(repository, reference.get("sdf")) if isinstance(reference, dict) else Path()
    if not sdf.is_file(): missing.append("Reference ligand SDF")
    mol2 = _profile_path(repository, reference.get("mol2")) if isinstance(reference, dict) else Path()
    if not mol2.is_file(): missing.append("Reference ligand MOL2")
    mapping = _profile_path(repository, reference.get("mapping")) if isinstance(reference, dict) else Path()
    if not mapping.is_file(): missing.append("Reference ligand chemistry mapping")
    for axis in "xyz":
        try: float(profile[f"center_{axis}"])
        except (KeyError, TypeError, ValueError): missing.append("Docking-box center"); break
    for axis in "xyz":
        try:
            if float(profile[f"size_{axis}"]) <= 0: raise ValueError
        except (KeyError, TypeError, ValueError): missing.append("Docking-box dimensions"); break
    try: find_vina_executable(repository.root)
    except FileNotFoundError: missing.append("Vina executable")
    return missing


def validate_ligand_pdbqt(path: Path) -> None:
    valid, reason = validate_pdbqt(path)
    if not valid: raise ValueError(reason)
    text = path.read_text(encoding="utf-8", errors="replace")
    for marker in ("ROOT", "TORSDOF"):
        if marker not in text: raise ValueError(f"Prepared ligand PDBQT is missing {marker}")
    atom_types = [line.split()[-1] for line in text.splitlines() if line.startswith(("ATOM", "HETATM"))]
    if not atom_types or any(not value.replace("+", "").replace("-", "").isalnum() for value in atom_types):
        raise ValueError("Prepared ligand PDBQT has invalid AutoDock atom types")


def validate_mol2(path: Path) -> None:
    if not path.is_file(): raise ValueError(f"Missing MOL2: {path.name}")
    text = path.read_text(encoding="utf-8", errors="replace")
    if "@<TRIPOS>ATOM" not in text or "@<TRIPOS>BOND" not in text: raise ValueError(f"Invalid MOL2 sections: {path.name}")
    header = text.splitlines(); marker = header.index("@<TRIPOS>MOLECULE")
    counts = header[marker + 2].split()
    if len(counts) < 2 or int(counts[0]) <= 0 or int(counts[1]) <= 0: raise ValueError(f"MOL2 contains no atoms or bonds: {path.name}")


class RedockingRunner:
    def __init__(self, repository: ProjectRepository, profile: dict[str, Any], settings: RedockingSettings,
                 *, progress: Callable[[dict[str, Any]], None] | None = None, cancel_event: threading.Event | None = None,
                 meeko: MeekoService | None = None, postdock: PostDockService | None = None,
                 vina: VinaDockingBackend | None = None) -> None:
        self.repository, self.profile, self.settings = repository, profile, settings
        self.progress = progress or (lambda event: None); self.cancel_event = cancel_event or threading.Event()
        self.meeko, self.postdock, self.vina = meeko or MeekoService(), postdock or PostDockService(), vina or VinaDockingBackend()

    def cancel(self) -> None: self.cancel_event.set()

    def _emit(self, event: str, stage: str, **values: Any) -> None: self.progress({"event": event, "stage": stage, **values})

    def _cancelled(self) -> bool: return self.cancel_event.is_set()

    def _stage_reusable(self, run_id: str, stage: str, fingerprint_value: str, paths: tuple[Path, ...]) -> bool:
        with self.repository.connection() as conn:
            row = conn.execute("SELECT status,fingerprint FROM redocking_stages WHERE run_id=? AND stage_name=?", (run_id, stage)).fetchone()
        if not (row and row["status"] == "completed" and row["fingerprint"] == fingerprint_value and all(path.is_file() and path.stat().st_size for path in paths)):
            return False
        try:
            if stage == "reference_validation": read_sdf(paths[0]); validate_mol2(paths[1])
            elif stage == "ligand_preparation": read_sdf(paths[0]); validate_ligand_pdbqt(paths[1])
            elif stage == "vina_redocking" and not parse_vina_poses(paths[0]): return False
            elif stage == "pose_export":
                for path in paths: read_sdf(path)
            elif stage == "mol2_generation":
                for path in paths: validate_mol2(path)
            elif stage == "final_validation":
                if json.loads(paths[0].read_text(encoding="utf-8")).get("status") != RedockingStatus.ARTIFACTS_READY.value: return False
        except Exception:
            return False
        return True

    def _stage_start(self, run_id: str, stage: str, fingerprint_value: str) -> None:
        if self._cancelled(): raise InterruptedError("Redocking cancelled")
        with self.repository.connection() as conn:
            conn.execute("UPDATE redocking_runs SET status=?,current_stage=?,error_stage=NULL,error_message=NULL WHERE run_id=?",
                         (RedockingStatus.RUNNING.value, stage, run_id))
            conn.execute("""INSERT INTO redocking_stages(run_id,stage_name,fingerprint,status,started_at)
                VALUES (?, ?, ?, 'running', ?) ON CONFLICT(run_id,stage_name) DO UPDATE SET
                fingerprint=excluded.fingerprint,status='running',started_at=excluded.started_at,finished_at=NULL,error_message=NULL""",
                (run_id, stage, fingerprint_value, utc_now()))
        self._emit("stage_started", stage)

    def _stage_done(self, run_id: str, stage: str, artifacts: dict[str, str]) -> None:
        with self.repository.connection() as conn:
            conn.execute("UPDATE redocking_stages SET status='completed',artifacts_json=?,finished_at=? WHERE run_id=? AND stage_name=?",
                         (json.dumps(artifacts, sort_keys=True), utc_now(), run_id, stage))
        self._emit("stage_completed", stage, artifacts=artifacts)

    def run(self, run_id: str | None = None) -> dict[str, Any]:
        missing = validate_redocking_prerequisites(self.repository, self.profile)
        if missing: raise ValueError("Redocking cannot start. Missing:\n• " + "\n• ".join(missing) + "\n\nReturn to Receptor Preparation to complete the profile.")
        profile_id = str(self.profile["id"]); reference_info = dict(self.profile["reference_ligand"])
        receptor = _profile_path(self.repository, self.profile["receptor"]); reference_sdf = _profile_path(self.repository, reference_info["sdf"])
        reference_mol2 = _profile_path(self.repository, reference_info["mol2"])
        executable = find_vina_executable(self.repository.root); run_id = run_id or uuid.uuid4().hex
        run_root = self.repository.root / "inputs" / "receptors" / profile_id / "redocking" / run_id
        input_dir, docking_dir, poses_dir, rmsd_dir = (run_root / name for name in ("input", "docking", "poses", "dockrmsd"))
        for folder in (input_dir, docking_dir, poses_dir, rmsd_dir): folder.mkdir(parents=True, exist_ok=True)
        settings = self.settings.as_dict(); receptor_hash, reference_hash = file_sha256(receptor), file_sha256(reference_sdf)
        with self.repository.connection() as conn:
            existing = conn.execute("SELECT run_id FROM redocking_runs WHERE run_id=?", (run_id,)).fetchone()
            if not existing:
                conn.execute("""INSERT INTO redocking_runs(run_id,receptor_profile_id,status,reference_ligand_id,receptor_path,receptor_sha256,
                    reference_sdf_path,reference_mol2_path,reference_sha256,settings_json,started_at)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", (run_id, profile_id, RedockingStatus.RUNNING.value,
                    str(reference_info["identity"]), receptor.relative_to(self.repository.root).as_posix(), receptor_hash,
                    reference_sdf.relative_to(self.repository.root).as_posix(), reference_mol2.relative_to(self.repository.root).as_posix(),
                    reference_hash, json.dumps(settings, sort_keys=True), utc_now()))
            else:
                conn.execute("""UPDATE redocking_runs SET status=?,interrupted_at=NULL,finished_at=NULL,cancel_requested=0,
                    receptor_path=?,receptor_sha256=?,reference_sdf_path=?,reference_mol2_path=?,reference_sha256=?,settings_json=? WHERE run_id=?""",
                    (RedockingStatus.RUNNING.value, receptor.relative_to(self.repository.root).as_posix(), receptor_hash,
                     reference_sdf.relative_to(self.repository.root).as_posix(), reference_mol2.relative_to(self.repository.root).as_posix(),
                     reference_hash, json.dumps(settings, sort_keys=True), run_id))
        try:
            ref_fp = fingerprint(settings={"mapping": "stable_atom_names"}, inputs={"sdf": reference_hash}, tool_version="reference-v1")
            if not self._stage_reusable(run_id, STAGES[0], ref_fp, (reference_sdf, reference_mol2)):
                self._stage_start(run_id, STAGES[0], ref_fp); reference = read_sdf(reference_sdf); validate_mol2(reference_mol2)
                mapping_path = _profile_path(self.repository, reference_info["mapping"]); mapping = json.loads(mapping_path.read_text(encoding="utf-8"))
                if int(mapping.get("heavy_atom_count", -1)) != int(chemistry_summary(reference)["heavy_atom_count"]): raise ValueError("Reference ligand atom mapping is incomplete")
                self._stage_done(run_id, STAGES[0], {"sdf": file_sha256(reference_sdf), "mol2": file_sha256(reference_mol2)})
            else: self._emit("stage_reused", STAGES[0])

            ligand_h, ligand_pdbqt = input_dir / "redocking_ligand_H.sdf", input_dir / "redocking_ligand.pdbqt"
            prep_fp = fingerprint(settings={"hydrogens": "rdkit_addHs", "arguments": ["-i", "-o"]}, inputs={"reference": reference_hash}, tool_version=_tool_version("meeko"))
            if not self._stage_reusable(run_id, STAGES[1], prep_fp, (ligand_h, ligand_pdbqt)):
                self._stage_start(run_id, STAGES[1], prep_fp); hydrogenated_copy(reference_sdf, ligand_h)
                result = self.meeko.prepare_ligand(ligand_h, ligand_pdbqt, run_root / "run.log"); validate_ligand_pdbqt(ligand_pdbqt)
                (input_dir / "meeko_command.json").write_text(json.dumps({"command": list(result.command), "exit_code": result.return_code,
                    "stdout": result.stdout, "stderr": result.stderr, "version": _tool_version("meeko")}, indent=2), encoding="utf-8")
                self._stage_done(run_id, STAGES[1], {"sdf": file_sha256(ligand_h), "pdbqt": file_sha256(ligand_pdbqt)})
            else: self._emit("stage_reused", STAGES[1])

            multi, vina_log = docking_dir / "redocked_multi_pose.pdbqt", docking_dir / "vina.log"
            vina_settings = {"center_x": self.profile["center_x"], "center_y": self.profile["center_y"], "center_z": self.profile["center_z"],
                "size_x": self.profile["size_x"], "size_y": self.profile["size_y"], "size_z": self.profile["size_z"],
                "exhaustiveness": settings["exhaustiveness"], "num_modes": settings["num_modes"], "energy_range": settings["energy_range"],
                "seed": settings["seed"], "cpu": settings["cpu_count"]}
            dock_fp = fingerprint(settings=vina_settings, inputs={"ligand": file_sha256(ligand_pdbqt), "receptor": receptor_hash}, tool_version=file_sha256(executable))
            if not self._stage_reusable(run_id, STAGES[2], dock_fp, (multi, vina_log)):
                self._stage_start(run_id, STAGES[2], dock_fp)
                result = self.vina.run(executable=executable, receptor=receptor, ligand=ligand_pdbqt, output=multi, log=vina_log,
                                       settings=vina_settings, cancelled=self._cancelled)
                if result.return_code or not result.poses: raise RuntimeError(result.stderr or "Vina produced no scored poses")
                (docking_dir / "command.json").write_text(json.dumps({"command": list(result.command), "exit_code": result.return_code}, indent=2), encoding="utf-8")
                self._stage_done(run_id, STAGES[2], {"pdbqt": file_sha256(multi), "log": file_sha256(vina_log)})
            else: self._emit("stage_reused", STAGES[2])
            scored = parse_vina_poses(multi)
            if not scored: raise ValueError("Vina output contains no scored poses")

            combined = poses_dir / "all_poses.sdf"
            export_fp = fingerprint(settings={"all_poses": True}, inputs={"pdbqt": file_sha256(multi)}, tool_version=_tool_version("meeko"))
            expected_sdfs = tuple(poses_dir / f"pose_{rank:03d}.sdf" for rank, *_ in scored)
            if not self._stage_reusable(run_id, STAGES[3], export_fp, expected_sdfs):
                self._stage_start(run_id, STAGES[3], export_fp); self.postdock.export_sdf(multi, combined, poses_dir / "meeko_export.log")
                Chem = __import__("rdkit.Chem", fromlist=["Chem"]); supplier = Chem.SDMolSupplier(str(combined), removeHs=False)
                molecules = [molecule for molecule in supplier if molecule is not None]
                if len(molecules) != len(scored): raise ValueError(f"Meeko exported {len(molecules)} structures for {len(scored)} scored poses")
                for (rank, affinity, _, _), molecule, destination in zip(scored, molecules, expected_sdfs):
                    molecule.SetIntProp("POSE_RANK", rank); molecule.SetDoubleProp("VINA_AFFINITY", affinity)
                    molecule.SetProp("RECEPTOR_PROFILE_ID", profile_id); molecule.SetProp("REDOCKING_RUN_ID", run_id)
                    molecule.SetIntProp("SEED", int(settings["seed"])); molecule.SetProp("BOX_JSON", json.dumps(vina_settings, sort_keys=True))
                    molecule.SetProp("SOURCE_PDBQT", multi.relative_to(self.repository.root).as_posix()); write_sdf(molecule, destination)
                self._stage_done(run_id, STAGES[3], {path.name: file_sha256(path) for path in expected_sdfs})
            else: self._emit("stage_reused", STAGES[3])

            reference_out = rmsd_dir / "reference_ligand_heavy.mol2"; mol2_paths = tuple(rmsd_dir / f"pose_{rank:03d}_heavy.mol2" for rank, *_ in scored)
            mol2_fp = fingerprint(settings={"hydrogens": "remove", "writer": "moldockpipe-mol2-v1"},
                                  inputs={path.name: file_sha256(path) for path in expected_sdfs}, tool_version="mol2-v1")
            if not self._stage_reusable(run_id, STAGES[4], mol2_fp, (reference_out, *mol2_paths)):
                self._stage_start(run_id, STAGES[4], mol2_fp); reference = read_sdf(reference_sdf)
                write_mol2(reference, reference_out, heavy_only=True, name="reference_ligand")
                for (rank, *_), sdf, mol2 in zip(scored, expected_sdfs, mol2_paths):
                    pose = read_sdf(sdf); validate_same_heavy_graph(reference, pose); write_mol2(pose, mol2, heavy_only=True, name=f"pose_{rank:03d}")
                shutil.copy2(mol2_paths[0], rmsd_dir / "top_ranked_pose_heavy.mol2")
                (rmsd_dir / "README.txt").write_text("DockRMSD is not run by MolDockPipe. Example commands:\n\n"
                    "DockRMSD reference_ligand_heavy.mol2 pose_001_heavy.mol2\n"
                    "DockRMSD reference_ligand_heavy.mol2 top_ranked_pose_heavy.mol2\n", encoding="utf-8")
                self._stage_done(run_id, STAGES[4], {path.name: file_sha256(path) for path in (reference_out, *mol2_paths)})
            else: self._emit("stage_reused", STAGES[4])

            final_fp = fingerprint(settings={"validation": "redocking-v1"}, inputs={"multi": file_sha256(multi), "reference": file_sha256(reference_out)}, tool_version="validation-v1")
            validation_path = run_root / "validation.json"
            if not self._stage_reusable(run_id, STAGES[5], final_fp, (validation_path,)):
                self._stage_start(run_id, STAGES[5], final_fp); validate_ligand_pdbqt(ligand_pdbqt); validate_mol2(reference_out)
                for path in (*expected_sdfs, *mol2_paths):
                    if not path.is_file() or not path.stat().st_size: raise ValueError(f"Missing pose artifact: {path.name}")
                for path in mol2_paths: validate_mol2(path)
                if file_sha256(rmsd_dir / "top_ranked_pose_heavy.mol2") != file_sha256(mol2_paths[0]): raise ValueError("Top-ranked MOL2 does not point to pose rank 1")
                pose_rows = []
                with self.repository.connection() as conn:
                    conn.execute("DELETE FROM redocking_poses WHERE run_id=?", (run_id,))
                    for (rank, affinity, _, _), sdf, mol2 in zip(scored, expected_sdfs, mol2_paths):
                        digest = file_sha256(mol2); metadata = {"seed": settings["seed"], "box": vina_settings, "source_pdbqt": multi.relative_to(self.repository.root).as_posix()}
                        conn.execute("INSERT INTO redocking_poses VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", (uuid.uuid4().hex, run_id, rank, affinity,
                            sdf.relative_to(self.repository.root).as_posix(), mol2.relative_to(self.repository.root).as_posix(),
                            multi.relative_to(self.repository.root).as_posix(), digest, json.dumps(metadata, sort_keys=True)))
                        pose_rows.append((rank, affinity, sdf.relative_to(self.repository.root).as_posix(), digest))
                with (run_root / "pose_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.writer(handle); writer.writerow(("pose_rank", "affinity", "sdf_path", "artifact_sha256")); writer.writerows(pose_rows)
                validation_path.write_text(json.dumps({"status": RedockingStatus.ARTIFACTS_READY.value, "rmsd_status": "Not calculated",
                    "pose_count": len(scored), "top_affinity": scored[0][1], "top_rank": 1, "hashes_recorded": True}, indent=2), encoding="utf-8")
                for path, artifact_type in ((ligand_h, "redocking_input_sdf"), (ligand_pdbqt, "redocking_ligand_pdbqt"),
                    (multi, "redocking_vina_output"), (vina_log, "redocking_vina_log"), (reference_out, "dockrmsd_reference_mol2"),
                    (rmsd_dir / "top_ranked_pose_heavy.mol2", "dockrmsd_top_pose_mol2"), (validation_path, "redocking_validation")):
                    self.repository.add_artifact(path, artifact_type, "redocking")
                for path in (*expected_sdfs, *mol2_paths):
                    self.repository.add_artifact(path, "redocking_pose_sdf" if path.suffix == ".sdf" else "dockrmsd_pose_mol2", "redocking")
                self._stage_done(run_id, STAGES[5], {"validation": file_sha256(validation_path)})
            else: self._emit("stage_reused", STAGES[5])
            with self.repository.connection() as conn:
                conn.execute("UPDATE redocking_runs SET status=?,finished_at=?,current_stage=NULL,prepared_ligand_path=?,prepared_ligand_sha256=?,meeko_version=?,vina_version=? WHERE run_id=?",
                    (RedockingStatus.ARTIFACTS_READY.value, utc_now(), ligand_pdbqt.relative_to(self.repository.root).as_posix(), file_sha256(ligand_pdbqt),
                     _tool_version("meeko"), executable.name, run_id))
            return {"run_id": run_id, "status": RedockingStatus.ARTIFACTS_READY.value, "run_root": str(run_root), "poses": len(scored), "top_affinity": scored[0][1]}
        except InterruptedError as exc:
            status = RedockingStatus.CANCELLED.value if self._cancelled() else RedockingStatus.INTERRUPTED.value
            with self.repository.connection() as conn:
                conn.execute("UPDATE redocking_runs SET status=?,interrupted_at=?,error_stage=current_stage,error_message=? WHERE run_id=?", (status, utc_now(), str(exc), run_id))
                conn.execute("UPDATE redocking_stages SET status=?,finished_at=?,error_message=? WHERE run_id=? AND status='running'", (status, utc_now(), str(exc), run_id))
            raise
        except Exception as exc:
            with self.repository.connection() as conn:
                row = conn.execute("SELECT current_stage FROM redocking_runs WHERE run_id=?", (run_id,)).fetchone(); stage = row[0] if row else None
                conn.execute("UPDATE redocking_runs SET status=?,finished_at=?,error_stage=?,error_message=? WHERE run_id=?", (RedockingStatus.FAILED.value, utc_now(), stage, str(exc)[:2000], run_id))
                conn.execute("UPDATE redocking_stages SET status='failed',finished_at=?,error_message=? WHERE run_id=? AND status='running'", (utc_now(), str(exc)[:2000], run_id))
            raise RuntimeError(f"Redocking failed during {stage}: {exc}") from exc
