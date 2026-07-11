from __future__ import annotations

import subprocess
import sys
from pathlib import Path


class PostDockService:
    def split(self, executable: Path, input_pdbqt: Path, output_dir: Path, prefix: str) -> list[Path]:
        output_dir.mkdir(parents=True, exist_ok=True)
        result = subprocess.run([str(executable), "--input", str(input_pdbqt), "--ligand", str(output_dir / prefix)], text=True, capture_output=True, check=False)
        if result.returncode:
            raise RuntimeError(result.stderr or result.stdout or "vina_split failed")
        return sorted(output_dir.glob(f"{prefix}*.pdbqt"))

    def export_sdf(self, input_pdbqt: Path, output_sdf: Path, log: Path) -> None:
        output_sdf.parent.mkdir(parents=True, exist_ok=True)
        log.parent.mkdir(parents=True, exist_ok=True)
        result = subprocess.run([sys.executable, "-m", "meeko.cli.mk_export", str(input_pdbqt), "-s", str(output_sdf)], text=True, capture_output=True, check=False)
        log.write_text(result.stdout + "\n" + result.stderr, encoding="utf-8")
        if result.returncode or not output_sdf.exists() or output_sdf.stat().st_size == 0:
            raise RuntimeError("mk_export.py failed")
