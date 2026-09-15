from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from collections.abc import Callable

from .executables import VINA, VINA_SPLIT, find_executable
from ..receptors.ad4zn import protocol_for, atom_types
from .validation import parse_vina_poses, validate_pdbqt


@dataclass(frozen=True)
class VinaResult:
    return_code: int
    poses: list[tuple[int, float, float, float]]
    stdout: str
    stderr: str
    command: tuple[str, ...] = ()


def find_vina_executable(project_root: Path) -> Path:
    return find_executable(VINA, project_root)


def find_vina_split_executable(project_root: Path) -> Path:
    return find_executable(VINA_SPLIT, project_root)


def build_vina_command(*, executable: Path, receptor: Path, ligand: Path, output: Path,
                       settings: dict[str, object], maps: Path | None = None) -> list[str]:
    protocol = protocol_for(settings)
    if protocol == "ad4zn":
        if maps is None:
            raise ValueError("AD4Zn docking requires validated AutoGrid maps")
        command = [str(executable), "--maps", str(maps), "--scoring", "ad4"]
        options = ("exhaustiveness", "num_modes", "energy_range", "seed", "cpu")
    else:
        command = [str(executable), "--receptor", str(receptor)]
        options = ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z",
                   "exhaustiveness", "num_modes", "energy_range", "seed", "cpu")
    command.extend(["--ligand", str(ligand), "--out", str(output)])
    for name in options:
        if name in settings and settings[name] is not None:
            command.extend([f"--{name}", str(settings[name])])
    return command


class VinaDockingBackend:
    def run(self, *, executable: Path, receptor: Path, ligand: Path, output: Path, log: Path, settings: dict[str, object],
            cancelled: Callable[[], bool] | None = None, maps: Path | None = None) -> VinaResult:
        valid, reason = validate_pdbqt(ligand)
        if not valid:
            raise ValueError(f"Ligand validation failed: {reason}")
        valid, reason = validate_pdbqt(receptor)
        if not valid:
            raise ValueError(f"Receptor validation failed: {reason}")
        output.parent.mkdir(parents=True, exist_ok=True)
        if protocol_for(settings) == "ad4zn":
            if maps is None:
                raise ValueError("AD4Zn maps are required")
            for suffix in (*atom_types(ligand), "e", "d"):
                path = Path(str(maps) + f".{suffix}.map")
                if not path.is_file() or not path.stat().st_size:
                    raise ValueError(f"Missing AD4Zn map for {suffix}")
        command = build_vina_command(executable=executable, receptor=receptor, ligand=ligand,
                                     output=output, settings=settings, maps=maps)
        process = subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=output.parent)
        while True:
            try:
                stdout, stderr = process.communicate(timeout=0.2)
                break
            except subprocess.TimeoutExpired:
                if cancelled and cancelled():
                    process.kill()
                    process.communicate()
                    raise InterruptedError("Vina docking was cancelled")
        log.write_text(
            "Vina command:\n"
            + " ".join(command)
            + "\n\nReturn code: "
            + str(process.returncode)
            + "\n\n--- stdout ---\n"
            + stdout
            + "\n--- stderr ---\n"
            + stderr,
            encoding="utf-8",
        )
        poses = parse_vina_poses(output) if process.returncode == 0 and output.exists() else []
        return VinaResult(int(process.returncode or 0), poses, stdout, stderr, tuple(command))
