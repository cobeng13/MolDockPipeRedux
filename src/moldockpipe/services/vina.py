from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

from .validation import parse_vina_poses, validate_pdbqt


@dataclass(frozen=True)
class VinaResult:
    return_code: int
    poses: list[tuple[int, float, float, float]]
    stdout: str
    stderr: str


class VinaDockingBackend:
    def run(self, *, executable: Path, receptor: Path, ligand: Path, output: Path, log: Path, settings: dict[str, object]) -> VinaResult:
        valid, reason = validate_pdbqt(ligand)
        if not valid:
            raise ValueError(f"Ligand validation failed: {reason}")
        valid, reason = validate_pdbqt(receptor)
        if not valid:
            raise ValueError(f"Receptor validation failed: {reason}")
        output.parent.mkdir(parents=True, exist_ok=True)
        command = [str(executable), "--receptor", str(receptor), "--ligand", str(ligand), "--out", str(output)]
        for name in ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z", "exhaustiveness", "num_modes", "energy_range", "seed", "cpu"):
            if name in settings and settings[name] is not None:
                command.extend([f"--{name}", str(settings[name])])
        completed = subprocess.run(command, text=True, capture_output=True, cwd=output.parent, check=False)
        log.write_text(
            "Vina command:\n"
            + " ".join(command)
            + "\n\nReturn code: "
            + str(completed.returncode)
            + "\n\n--- stdout ---\n"
            + completed.stdout
            + "\n--- stderr ---\n"
            + completed.stderr,
            encoding="utf-8",
        )
        poses = parse_vina_poses(output) if completed.returncode == 0 and output.exists() else []
        return VinaResult(completed.returncode, poses, completed.stdout, completed.stderr)
