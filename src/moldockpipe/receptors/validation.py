from __future__ import annotations

import math
from pathlib import Path

from .models import ReceptorPreparationPlan


def validate_plan(plan: ReceptorPreparationPlan) -> None:
    numbers = (*plan.box_center, *plan.box_size)
    if not all(math.isfinite(value) for value in numbers):
        raise ValueError("Search-box coordinates and dimensions must be finite")
    if any(value <= 0 for value in plan.box_size):
        raise ValueError("Search-box dimensions must be positive")
    if not plan.included_chains:
        raise ValueError("At least one receptor chain must be included")


def validate_prepared_outputs(folder: Path, plan: ReceptorPreparationPlan) -> list[str]:
    warnings: list[str] = []
    pdbqt = folder / "receptor.pdbqt"
    if not pdbqt.is_file() or pdbqt.stat().st_size == 0:
        raise RuntimeError("Meeko did not create a nonempty receptor PDBQT")
    text = pdbqt.read_text(encoding="utf-8", errors="replace")
    if not any(line.startswith(("ATOM", "HETATM")) for line in text.splitlines()):
        raise RuntimeError("Prepared receptor PDBQT contains no receptor atoms")
    if plan.reference_ligand:
        reference = folder / "reference_ligand" / "reference_ligand.pdb"
        coordinates = []
        for line in reference.read_text(encoding="ascii").splitlines():
            if line.startswith(("ATOM", "HETATM")):
                coordinates.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        for point in coordinates:
            if any(abs(point[i] - plan.box_center[i]) > plan.box_size[i] / 2 + 1e-6 for i in range(3)):
                warnings.append("Reference ligand extends outside the selected docking box")
                break
    return warnings
