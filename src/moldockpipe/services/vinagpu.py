from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from collections.abc import Callable

from ..fingerprints import file_sha256, fingerprint
from .executables import application_root
from .validation import parse_vina_poses, validate_pdbqt
from .vina import VinaResult


KERNEL_NAMES = ("Kernel1_Opt.bin", "Kernel2_Opt.bin")


def find_vinagpu_executable(project_root: Path) -> Path:
    """Find Vina-GPU in tools/vinagpu, unless explicitly overridden."""
    override = os.environ.get("MOLDOCKPIPE_VINAGPU", "").strip()
    names = ("Vina-GPU+.exe", "Vina-GPU.exe", "vina-gpu.exe", "vina-gpu")
    if override:
        candidate = Path(override).expanduser()
        located = str(candidate) if candidate.is_file() else shutil.which(override)
        if located:
            return Path(located).resolve()
        raise FileNotFoundError(f"MOLDOCKPIPE_VINAGPU does not point to a Vina-GPU executable: {override}")
    for root in (application_root(), project_root.resolve()):
        folder = root / "tools" / "vinagpu"
        for name in names:
            candidate = folder / name
            if candidate.is_file() and (os.name == "nt" or os.access(candidate, os.X_OK)):
                return candidate.resolve()
    raise FileNotFoundError("Vina-GPU executable not found. Set MOLDOCKPIPE_VINAGPU or install it with its kernels in tools/vinagpu.")


def vinagpu_tool_hash(executable: Path) -> str:
    kernels = {}
    for name in KERNEL_NAMES:
        path = executable.parent / name
        if not path.is_file():
            raise FileNotFoundError(f"Vina-GPU kernel missing beside executable: {path}")
        kernels[name] = file_sha256(path)
    return fingerprint(settings={"backend": "vinagpu", "kernels": kernels}, inputs={},
                       tool_version=file_sha256(executable))


class VinaGpuDockingBackend:
    """Single-state GPU docking with isolated kernel snapshots and normal Vina artifacts."""

    def run(self, *, executable: Path, receptor: Path, ligand: Path, output: Path, log: Path,
            settings: dict[str, object], cancelled: Callable[[], bool] | None = None) -> VinaResult:
        for label, path in (("Ligand", ligand), ("Receptor", receptor)):
            valid, reason = validate_pdbqt(path)
            if not valid:
                raise ValueError(f"{label} validation failed: {reason}")
        thread = int(settings.get("thread", 8000))
        depth = int(settings.get("search_depth", 32))
        if thread < 1000 or depth < 1:
            raise ValueError("Vina-GPU requires at least 1000 computing lanes and a positive search depth")
        executable = executable.resolve()
        vinagpu_tool_hash(executable)
        output = output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        for name in KERNEL_NAMES:
            shutil.copy2(executable.parent / name, output.parent / name)
        command = [str(executable), "--receptor", str(receptor.resolve()),
                   "--ligand", str(ligand.resolve()), "--out", str(output),
                   "--thread", str(thread), "--search_depth", str(depth)]
        for name in ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z", "num_modes", "energy_range", "seed"):
            if name in settings and settings[name] is not None:
                command.extend([f"--{name}", str(settings[name])])
        if cancelled and cancelled():
            raise InterruptedError("Vina-GPU docking was cancelled")
        process = subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   cwd=output.parent)
        interrupted = False
        while True:
            try:
                stdout, stderr = process.communicate(timeout=0.1)
                break
            except subprocess.TimeoutExpired:
                if cancelled and cancelled():
                    interrupted = True
                    process.terminate()
                    try:
                        stdout, stderr = process.communicate(timeout=5)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        stdout, stderr = process.communicate()
                    break
        log.write_text("Vina-GPU command:\n" + " ".join(command)
                       + f"\n\nReturn code: {process.returncode}\n\n--- stdout ---\n"
                       + stdout + "\n--- stderr ---\n" + stderr, encoding="utf-8")
        if interrupted:
            raise InterruptedError("Vina-GPU docking was cancelled")
        valid, _ = validate_pdbqt(output)
        poses = parse_vina_poses(output) if process.returncode == 0 and valid else []
        return VinaResult(int(process.returncode or 0), poses, stdout, stderr, tuple(command))
