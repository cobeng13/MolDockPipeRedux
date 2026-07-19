from __future__ import annotations

from pathlib import Path

from .models import ReceptorPreparationPlan, ResidueKey


def structure_as_pdb(path: Path) -> str:
    try:
        import gemmi
    except ImportError as exc:
        raise RuntimeError("Gemmi is required to read receptor structures") from exc
    structure = gemmi.read_structure(str(path))
    return structure.make_pdb_string()


def _line_key(line: str) -> ResidueKey:
    return ResidueKey(line[21:22].strip(), int(line[22:26]), line[26:27].strip(), line[17:20].strip())


def _selected_atom_lines(pdb_text: str, plan: ReceptorPreparationPlan) -> list[str]:
    candidates: dict[ResidueKey, list[str]] = {}
    current_model = 0
    explicit_models = False
    for line in pdb_text.splitlines():
        record = line[:6].strip()
        if record == "MODEL":
            explicit_models = True
            try:
                current_model = int(line[10:14].strip()) - 1
            except ValueError:
                current_model += 1
            continue
        if record == "ENDMDL":
            continue
        if record not in {"ATOM", "HETATM"}:
            continue
        if explicit_models and current_model != plan.selected_model:
            continue
        key = _line_key(line)
        if plan.included_chains and key.chain not in plan.included_chains:
            continue
        candidates.setdefault(key, []).append(line)
    selected: list[str] = []
    for key, lines in candidates.items():
        specified = plan.altloc_choices.get(key)
        if specified is None:
            # Prefer the conformer with the largest atom set.  mmCIF files can
            # legitimately contain only B/C for a residue, so hard-coding A
            # silently creates an incomplete amino acid.
            alternatives = sorted({line[16:17].strip() for line in lines if line[16:17].strip()})
            if alternatives:
                def score(alt: str) -> tuple[int, float, bool]:
                    active = [line for line in lines if not line[16:17].strip() or line[16:17].strip() == alt]
                    return (len({line[12:16].strip() for line in active}), sum(float(line[54:60] or 0) for line in active), alt == "A")
                selected_altloc = max(alternatives, key=score)
            else:
                selected_altloc = ""
        else:
            selected_altloc = specified
        for line in lines:
            altloc = line[16:17].strip()
            if altloc and altloc != selected_altloc:
                continue
            selected.append(line[:16] + " " + line[17:] if altloc else line)
    return selected


def write_cleaned_receptor(pdb_text: str, plan: ReceptorPreparationPlan, destination: Path) -> None:
    removed = set(plan.removed_residues) | set(plan.excluded_receptor_residues)
    if plan.reference_ligand:
        removed.add(plan.reference_ligand)
    lines = [line for line in _selected_atom_lines(pdb_text, plan) if _line_key(line) not in removed]
    if not lines:
        raise ValueError("The receptor composition removed every atom")
    destination.write_text("\n".join(lines + ["END", ""]), encoding="ascii")


def write_reference_ligand(pdb_text: str, plan: ReceptorPreparationPlan, destination: Path) -> None:
    if not plan.reference_ligand:
        return
    lines = [line for line in _selected_atom_lines(pdb_text, plan) if _line_key(line) == plan.reference_ligand]
    if not lines:
        raise ValueError(f"Reference ligand {plan.reference_ligand.label()} was not found")
    destination.write_text("\n".join(lines + ["END", ""]), encoding="ascii")


def write_box_pdb(center: tuple[float, float, float], size: tuple[float, float, float], destination: Path) -> None:
    corners = []
    for dx in (-0.5, 0.5):
        for dy in (-0.5, 0.5):
            for dz in (-0.5, 0.5):
                corners.append((center[0] + dx*size[0], center[1] + dy*size[1], center[2] + dz*size[2]))
    lines = [f"HETATM{i:5d}  C{i} BOX X   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
             for i, (x, y, z) in enumerate(corners, 1)]
    destination.write_text("\n".join(lines + ["END", ""]), encoding="ascii")
