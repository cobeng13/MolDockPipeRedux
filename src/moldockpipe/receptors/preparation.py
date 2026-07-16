from __future__ import annotations

import hashlib
import json
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


def prepare_receptor(project_root: Path, plan: ReceptorPreparationPlan) -> Path:
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
            raise RuntimeError(f"Meeko receptor preparation failed (exit {result.returncode}). See preparation log.\n{result.stderr[-1000:]}")
        warnings = validate_prepared_outputs(work, plan)
        report: dict[str, Any] = {
            "source_sha256": _sha256(source), "cleaned_receptor_sha256": _sha256(cleaned),
            "prepared_receptor_sha256": _sha256(work / "receptor.pdbqt"), "meeko_version": _meeko_version(),
            "created_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            "chains_included": list(plan.included_chains),
            "components_removed": [key.label() for key in plan.removed_residues],
            "components_retained": [key.label() for key in plan.retained_components],
            "reference_ligand": plan.reference_ligand.label() if plan.reference_ligand else None,
            "center": dict(zip(("x", "y", "z"), plan.box_center)),
            "box": {"method": plan.box_method, **plan.box_parameters,
                    **dict(zip(("size_x", "size_y", "size_z"), plan.box_size))},
            "command": command, "return_code": result.returncode, "warnings": warnings,
        }
        (work / "preparation_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        work.rename(final)
    return final
