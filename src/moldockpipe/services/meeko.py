from __future__ import annotations

import subprocess
import sys
from pathlib import Path


class MeekoService:
    """CLI adapter kept behind a service boundary for Meeko API compatibility."""

    def prepare_ligand(self, sdf: Path, pdbqt: Path, log: Path) -> None:
        pdbqt.parent.mkdir(parents=True, exist_ok=True)
        log.parent.mkdir(parents=True, exist_ok=True)
        command = [sys.executable, "-m", "meeko.cli.mk_prepare_ligand", "-i", str(sdf), "-o", str(pdbqt)]
        result = subprocess.run(command, text=True, capture_output=True, check=False)
        log.write_text(result.stdout + "\n" + result.stderr, encoding="utf-8")
        if result.returncode or not pdbqt.exists() or pdbqt.stat().st_size == 0:
            raise RuntimeError("Meeko ligand preparation failed")

    def export_docked_poses(self, vina_pdbqt: Path, output_sdf: Path, log: Path) -> None:
        # Meeko's export interface varies by release; probe the documented module form.
        command = [sys.executable, "-m", "meeko.cli_export", "-i", str(vina_pdbqt), "-o", str(output_sdf)]
        output_sdf.parent.mkdir(parents=True, exist_ok=True)
        log.parent.mkdir(parents=True, exist_ok=True)
        result = subprocess.run(command, text=True, capture_output=True, check=False)
        log.write_text(result.stdout + "\n" + result.stderr, encoding="utf-8")
        if result.returncode or not output_sdf.exists() or output_sdf.stat().st_size == 0:
            raise RuntimeError("Meeko docked-pose export failed; verify the installed Meeko export command")
