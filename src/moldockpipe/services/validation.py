from __future__ import annotations

import math
import re
from pathlib import Path


VINA_RESULT = re.compile(r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)")


def validate_pdbqt(path: Path) -> tuple[bool, str]:
    if not path.exists() or path.stat().st_size < 100:
        return False, "PDBQT is missing or too small"
    atoms = [line for line in path.read_text(errors="replace").splitlines() if line.startswith(("ATOM", "HETATM"))]
    if not atoms:
        return False, "PDBQT has no atom records"
    for atom in atoms:
        try:
            coordinates = [float(atom[30:38]), float(atom[38:46]), float(atom[46:54])]
        except ValueError:
            return False, "PDBQT contains malformed coordinates"
        if not all(math.isfinite(value) for value in coordinates):
            return False, "PDBQT contains non-finite coordinates"
    return True, "valid"


def parse_vina_poses(path: Path) -> list[tuple[int, float, float, float]]:
    poses = []
    for mode_index, match in enumerate(VINA_RESULT.finditer(path.read_text(errors="replace")), 1):
        poses.append((mode_index, float(match.group(1)), float(match.group(2)), float(match.group(3))))
    return poses
