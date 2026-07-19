from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

from .extraction import structure_as_pdb, write_box_pdb, write_cleaned_receptor, write_reference_ligand
from .models import ReceptorPreparationPlan
from .ligand_chemistry import create_reference_bundle
from .validation import validate_plan, validate_prepared_outputs
from .analysis import analysis_as_record, analyze_structure
from ..models import StageStatus
from ..project import ProjectRepository


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _meeko_version() -> str:
    try:
        from importlib.metadata import version
        return version("meeko")
    except Exception:
        return "unknown"


def _preserve_failure_artifacts(project_root: Path, profile_id: str, work: Path) -> Path:
    """Keep the exact Meeko input/log after a failed temporary preparation."""
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    destination = Path(project_root) / "logs" / "receptor_preparation" / profile_id / timestamp
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(work, destination)
    return destination


def _meeko_failure_message(stderr: str, failure_dir: Path) -> str:
    residues = sorted(set(re.findall(r"matched with excess inter-residue bond\(s\):\s*([^\s]+)", stderr)))
    if residues:
        summary = ", ".join(residues)
        reason = (
            f"Meeko detected unexpected inter-residue bonds near {summary}. "
            "This usually indicates a structural clash, unresolved alternate locations, or a covalent/modified residue."
        )
    else:
        reason = "Meeko could not match one or more receptor residues to its preparation templates."
    return (
        "Meeko receptor preparation failed.\n\n"
        f"{reason}\n\n"
        "The wizard did not delete residues automatically. Review the saved cleaned receptor and log, then either choose a different alternate location, "
        "remove the affected component deliberately, or repair the source structure externally.\n\n"
        f"Failure artifacts: {failure_dir}"
    )


def _prepare_receptor_files(project_root: Path, plan: ReceptorPreparationPlan,
                            inventory: dict[str, object]) -> Path:
    """Prepare into a temporary sibling and publish only after every check passes."""
    validate_plan(plan)
    parent = Path(project_root) / "inputs" / "receptors"
    parent.mkdir(parents=True, exist_ok=True)
    final = parent / plan.profile_id
    if final.exists():
        raise FileExistsError(f"Receptor profile folder already exists: {final}")
    with tempfile.TemporaryDirectory(prefix=f".{plan.profile_id}-", dir=parent) as temporary:
        work = Path(temporary)
        source_suffix = plan.source_path.suffix.lower() or ".pdb"
        source = work / f"source{source_suffix}"
        shutil.copy2(plan.source_path, source)
        pdb_text = structure_as_pdb(source)
        inventory_path = work / "structure_inventory.json"
        inventory_path.write_text(json.dumps(inventory, indent=2), encoding="utf-8")
        cleaned = work / "receptor_cleaned.pdb"
        write_cleaned_receptor(pdb_text, plan, cleaned)
        if plan.reference_ligand:
            reference_folder = work / "reference_ligand"
            reference_folder.mkdir()
            reference_pdb = reference_folder / "reference_ligand.pdb"
            write_reference_ligand(pdb_text, plan, reference_pdb)
            template = None
            if plan.chemistry_template_path:
                template = reference_folder / "chemical_template.sdf"
                shutil.copy2(plan.chemistry_template_path, template)
            create_reference_bundle(reference_pdb, reference_folder, template_sdf=template)
        write_box_pdb(plan.box_center, plan.box_size, work / "receptor.box.pdb")
        (work / "preparation.yml").write_text(yaml.safe_dump(plan.as_record(), sort_keys=False), encoding="utf-8")

        command = [sys.executable, "-m", "meeko.cli.mk_prepare_receptor", "--read_pdb", str(cleaned),
                   "--output_basename", "receptor", "--write_pdbqt", "receptor.pdbqt",
                   "--write_json", "receptor.json", "--write_vina_box", "receptor.box.txt",
                   "--write_pdb", "receptor_prepared.pdb", "--box_center", *(f"{v:.6f}" for v in plan.box_center),
                   "--box_size", *(f"{v:.6f}" for v in plan.box_size), "--default_altloc", "A"]
        result = subprocess.run(command, cwd=work, text=True, capture_output=True, check=False)
        log_text = "$ " + subprocess.list2cmdline(command) + "\n\n" + result.stdout + "\n" + result.stderr
        (work / "preparation.log").write_text(log_text, encoding="utf-8")
        if result.returncode:
            failure_dir = _preserve_failure_artifacts(project_root, plan.profile_id, work)
            raise RuntimeError(_meeko_failure_message(result.stderr, failure_dir))
        warnings = validate_prepared_outputs(work, plan)
        report: dict[str, Any] = {
            "receptor_profile_id": plan.profile_id, "receptor_name": plan.profile_name,
            "source_sha256": _sha256(source), "cleaned_receptor_sha256": _sha256(cleaned),
            "prepared_receptor_sha256": _sha256(work / "receptor.pdbqt"), "meeko_version": _meeko_version(),
            "structure_inventory_sha256": _sha256(inventory_path),
            "created_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            "chains_included": list(plan.included_chains),
            "components_removed": [key.label() for key in plan.removed_residues],
            "components_retained": [key.label() for key in plan.retained_components],
            "receptor_residues_excluded": [key.label() for key in plan.excluded_receptor_residues],
            "altloc_choices": {key.label(): value for key, value in plan.altloc_choices.items()},
            "reference_ligand": plan.reference_ligand.label() if plan.reference_ligand else None,
            "center": {"method": plan.center_method, "parameters": plan.center_parameters,
                       **dict(zip(("x", "y", "z"), plan.box_center))},
            "box": {"method": plan.box_method, **plan.box_parameters,
                    **dict(zip(("size_x", "size_y", "size_z"), plan.box_size))},
            "command": command, "return_code": result.returncode, "warnings": warnings,
        }
        (work / "preparation_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        work.rename(final)
    return final


def prepare_receptor(project_root: Path, plan: ReceptorPreparationPlan) -> Path:
    """Prepare a receptor and persist report-ready provenance for success or failure."""
    project_root = Path(project_root)
    repository = ProjectRepository(project_root)
    analysis = analyze_structure(plan.source_path, plan.selected_model, plan.included_chains)
    inventory = analysis_as_record(analysis, selected_model=plan.selected_model,
                                   included_chains=plan.included_chains)
    preparation_run_id = repository.start_receptor_preparation(
        profile_id=plan.profile_id,
        receptor_name=plan.profile_name,
        source_path=plan.source_path,
        source_sha256=_sha256(plan.source_path),
        inventory=inventory,
        decisions=plan.as_record(),
    )
    try:
        final = _prepare_receptor_files(project_root, plan, inventory)
    except Exception as exc:
        repository.finish_receptor_preparation(
            preparation_run_id, StageStatus.FAILED.value, error=str(exc)[:4000]
        )
        raise
    report = json.loads((final / "preparation_report.json").read_text(encoding="utf-8"))
    artifacts = {path.relative_to(final).as_posix(): {
                    "path": path.relative_to(project_root).as_posix(), "sha256": _sha256(path),
                 } for path in sorted(final.rglob("*")) if path.is_file()}
    repository.finish_receptor_preparation(
        preparation_run_id,
        StageStatus.COMPLETED.value,
        command={"argv": report.get("command", []), "return_code": report.get("return_code")},
        artifacts=artifacts,
        warnings=list(report.get("warnings", [])),
    )
    return final
