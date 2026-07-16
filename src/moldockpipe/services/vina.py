from __future__ import annotations

import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from collections.abc import Callable

from .validation import parse_vina_poses, validate_pdbqt


@dataclass(frozen=True)
class VinaResult:
    return_code: int
    poses: list[tuple[int, float, float, float]]
    stdout: str
    stderr: str
    command: tuple[str, ...] = ()


def find_vina_executable(project_root: Path) -> Path:
    application_root = Path(__file__).resolve().parents[3]
    candidates = (
        application_root / "tools" / "vina" / "vina.exe", application_root / "tools" / "vina" / "vina_1.2.7_win.exe",
        project_root / "tools" / "vina" / "vina.exe", project_root / "tools" / "vina" / "vina_1.2.7_win.exe",
    )
    executable = next((candidate for candidate in candidates if candidate.is_file()), None)
    if executable is None: raise FileNotFoundError(f"Vina executable not found. Place it at {application_root / 'tools' / 'vina'}")
    return executable


class VinaDockingBackend:
    def run(self, *, executable: Path, receptor: Path, ligand: Path, output: Path, log: Path, settings: dict[str, object],
            cancelled: Callable[[], bool] | None = None) -> VinaResult:
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
        process = subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=output.parent)
        while process.poll() is None:
            if cancelled and cancelled():
                process.terminate()
                try: process.wait(timeout=5)
                except subprocess.TimeoutExpired: process.kill()
                stdout, stderr = process.communicate()
                raise InterruptedError("Vina redocking was cancelled")
            time.sleep(0.05)
        stdout, stderr = process.communicate()
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
