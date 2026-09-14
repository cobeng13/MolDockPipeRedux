from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path


def receptor_export_filename(profile_name: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", profile_name.strip()).strip("._") or "receptor"
    return f"{safe}_prepared.pdb"


def _pdbqt_element(line: str) -> str:
    atom_type = line[76:].split()[-1] if line[76:].split() else ""
    mapped = {"A": "C", "NA": "N", "OA": "O", "SA": "S", "HD": "H", "HS": "H"}.get(atom_type)
    if mapped:
        return mapped
    letters = "".join(character for character in atom_type if character.isalpha())
    if letters:
        return letters[:2].capitalize()
    atom_name = line[12:16].strip()
    letters = "".join(character for character in atom_name if character.isalpha())
    return (letters[:2] if len(letters) > 1 and letters[1].islower() else letters[:1]).capitalize()


class PostDockService:
    def export_validation_standard(self, source_mol2: Path, output_mol2: Path) -> None:
        """Copy the exact heavy-atom MOL2 standard from a completed validation."""
        if not source_mol2.is_file():
            raise FileNotFoundError(f"Validation reference MOL2 not found: {source_mol2}")
        text = source_mol2.read_text(encoding="utf-8", errors="replace")
        for section in ("@<TRIPOS>MOLECULE", "@<TRIPOS>ATOM", "@<TRIPOS>BOND"):
            if section not in text:
                raise ValueError(f"Validation reference MOL2 is missing {section}: {source_mol2}")
        output_mol2.parent.mkdir(parents=True, exist_ok=True)
        temporary = output_mol2.with_suffix(output_mol2.suffix + ".tmp")
        shutil.copy2(source_mol2, temporary)
        temporary.replace(output_mol2)

    def export_receptor_pdb(self, receptor_pdbqt: Path, output_pdb: Path,
                            prepared_pdb: Path | None = None) -> str:
        """Export the post-Meeko receptor PDB, with a PDBQT-derived fallback.

        The cleaned receptor is deliberately not used: it is the input to
        Meeko and may not contain the protonation/atom representation used for
        docking.
        """
        output_pdb.parent.mkdir(parents=True, exist_ok=True)
        temporary = output_pdb.with_suffix(output_pdb.suffix + ".tmp")
        if prepared_pdb and prepared_pdb.is_file():
            text = prepared_pdb.read_text(encoding="utf-8", errors="replace")
            if any(line.startswith(("ATOM  ", "HETATM")) for line in text.splitlines()):
                shutil.copy2(prepared_pdb, temporary)
                temporary.replace(output_pdb)
                return "meeko_prepared_pdb"
        if not receptor_pdbqt.is_file():
            raise FileNotFoundError(f"Prepared receptor PDBQT not found: {receptor_pdbqt}")
        atoms = []
        for line in receptor_pdbqt.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                atoms.append(line[:66].ljust(76) + _pdbqt_element(line).rjust(2))
            elif line.startswith("TER"):
                atoms.append("TER")
        if not atoms:
            raise ValueError(f"Prepared receptor PDBQT contains no receptor atoms: {receptor_pdbqt}")
        temporary.write_text(
            "REMARK 950 Generated from the prepared receptor PDBQT used for docking\n"
            + "\n".join(atoms) + "\nEND\n", encoding="ascii"
        )
        temporary.replace(output_pdb)
        return "prepared_pdbqt_conversion"

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
