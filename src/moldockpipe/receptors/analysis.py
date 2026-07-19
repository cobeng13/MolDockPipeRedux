from __future__ import annotations

from pathlib import Path

from .classification import WATER_NAMES, classify_component
from .models import AtomRecord, Component, ProteinResidueIssue, ResidueKey, StructureAnalysis


# Heavy atoms expected for the conventional amino-acid residues Meeko handles.
# This is deliberately limited to standard residues: modified residues remain a
# Meeko/template concern rather than being guessed by the wizard.
_STANDARD_RESIDUE_ATOMS = {
    "ALA": "N CA C O CB", "ARG": "N CA C O CB CG CD NE CZ NH1 NH2",
    "ASN": "N CA C O CB CG OD1 ND2", "ASP": "N CA C O CB CG OD1 OD2",
    "CYS": "N CA C O CB SG", "GLN": "N CA C O CB CG CD OE1 NE2",
    "GLU": "N CA C O CB CG CD OE1 OE2", "GLY": "N CA C O",
    "HIS": "N CA C O CB CG ND1 CD2 CE1 NE2", "ILE": "N CA C O CB CG1 CG2 CD1",
    "LEU": "N CA C O CB CG CD1 CD2", "LYS": "N CA C O CB CG CD CE NZ",
    "MET": "N CA C O CB CG SD CE", "PHE": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ",
    "PRO": "N CA C O CB CG CD", "SER": "N CA C O CB OG",
    "THR": "N CA C O CB OG1 CG2", "TRP": "N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
    "TYR": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH", "VAL": "N CA C O CB CG1 CG2",
}


def _best_altloc(atoms: tuple[AtomRecord, ...]) -> str:
    """Choose the most complete conformer, preferring A only on a tie."""
    options = sorted({atom.altloc for atom in atoms if atom.altloc})
    if not options:
        return ""
    unlabelled = {atom.name for atom in atoms if not atom.altloc}
    def score(label: str) -> tuple[int, float, int]:
        chosen = [atom for atom in atoms if not atom.altloc or atom.altloc == label]
        return (len(unlabelled | {atom.name for atom in chosen}), sum(atom.occupancy for atom in chosen), label == "A")
    return max(options, key=score)


def _protein_issue(key: ResidueKey, atoms: tuple[AtomRecord, ...]) -> ProteinResidueIssue | None:
    expected = _STANDARD_RESIDUE_ATOMS.get(key.name)
    if not expected:
        return None
    altloc = _best_altloc(atoms)
    selected = {atom.name for atom in atoms if not atom.altloc or atom.altloc == altloc}
    missing = tuple(sorted(set(expected.split()) - selected))
    if not missing:
        return None
    locations = tuple(sorted({atom.altloc for atom in atoms if atom.altloc}))
    return ProteinResidueIssue(key, missing, locations, altloc)


def analysis_as_record(analysis: StructureAnalysis, *, selected_model: int,
                       included_chains: tuple[str, ...]) -> dict[str, object]:
    """Create a compact, report-ready snapshot of the unmodified structure."""
    return {
        "source_path": str(analysis.source_path),
        "models": list(analysis.models),
        "selected_model": selected_model,
        "chains_by_model": {str(key): list(value) for key, value in analysis.chains_by_model.items()},
        "included_chains": list(included_chains),
        "counts": {
            "protein_residues": analysis.protein_residue_count,
            "waters": analysis.water_count,
            "nonpolymer_components": len(analysis.components),
            "alternate_location_atoms": len(analysis.alternate_locations),
            "zero_occupancy_atoms": len(analysis.zero_occupancy_atoms),
            "insertion_codes": len(analysis.insertion_codes),
        },
        "protein_residues": [key.label() for key in analysis.protein_residues],
        "components": [{
            "identity": component.key.label(),
            "name": component.key.name,
            "chain": component.key.chain,
            "residue_number": component.key.number,
            "insertion_code": component.key.insertion_code,
            "atom_count": len(component.atoms),
            "heavy_atom_count": len(component.heavy_atoms),
            "category": component.category,
            "suggested_role": component.suggested_role.value,
            "classification_reason": component.reason,
            "covalently_linked": component.covalently_linked,
        } for component in analysis.components],
        "alternate_locations": list(analysis.alternate_locations),
        "zero_occupancy_atoms": list(analysis.zero_occupancy_atoms),
        "insertion_code_residues": list(analysis.insertion_codes),
        "connections": list(analysis.connections),
        "incomplete_protein_residues": [{
            "identity": issue.key.label(),
            "missing_atoms": list(issue.missing_atoms),
            "alternate_locations": list(issue.alternate_locations),
            "recommended_altloc": issue.recommended_altloc,
        } for issue in analysis.protein_residue_issues],
    }


def _gemmi():
    try:
        import gemmi
    except ImportError as exc:
        raise RuntimeError("Gemmi is required to analyze receptor structures. Install the project dependencies first.") from exc
    return gemmi


def analyze_structure(path: Path, model_index: int = 0, included_chains: tuple[str, ...] = ()) -> StructureAnalysis:
    path = Path(path)
    if path.suffix.lower() not in {".pdb", ".cif", ".mmcif"}:
        raise ValueError("Structure must be a .pdb, .cif, or .mmcif file")
    gemmi = _gemmi()
    structure = gemmi.read_structure(str(path))
    if len(structure) == 0:
        raise ValueError("Structure contains no models")
    if model_index < 0 or model_index >= len(structure):
        raise ValueError(f"Model index {model_index + 1} is not present")
    models = tuple(str(getattr(model, "name", index + 1)) for index, model in enumerate(structure))
    chains_by_model = {index: tuple(chain.name for chain in model) for index, model in enumerate(structure)}
    selected = structure[model_index]
    allowed = set(included_chains)
    components: list[Component] = []
    protein_count = water_count = 0
    protein_issues: list[ProteinResidueIssue] = []
    protein_residues: list[ResidueKey] = []
    altlocs: set[str] = set()
    insertion_codes: set[str] = set()
    zero_occupancy: list[str] = []

    for chain in selected:
        if allowed and chain.name not in allowed:
            continue
        for residue in chain:
            name = residue.name.strip().upper()
            seqid = residue.seqid
            icode = str(seqid.icode).strip()
            key = ResidueKey(chain.name, int(seqid.num), icode, name)
            atoms = tuple(AtomRecord(
                atom.name.strip(), atom.element.name, float(atom.pos.x), float(atom.pos.y), float(atom.pos.z),
                float(atom.occ), str(atom.altloc).strip().replace("\x00", ""),
            ) for atom in residue)
            for atom in atoms:
                if atom.altloc:
                    altlocs.add(f"{key.label()}:{atom.name}:{atom.altloc}")
                if atom.occupancy <= 0:
                    zero_occupancy.append(f"{key.label()}:{atom.name}")
            if icode:
                insertion_codes.add(key.label())
            info = gemmi.find_tabulated_residue(name)
            amino = bool(info and info.is_amino_acid())
            water = name in WATER_NAMES or bool(info and info.is_water())
            # Gemmi uses 'A' for ATOM and 'H' for HETATM records.
            if amino and str(getattr(residue, "het_flag", "A")) != "H":
                protein_count += 1
                protein_residues.append(key)
                issue = _protein_issue(key, atoms)
                if issue:
                    protein_issues.append(issue)
                continue
            if water:
                water_count += 1
            polymer_like = amino
            category, role, reason = classify_component(name, atoms, polymer_like=polymer_like)
            components.append(Component(key, atoms, role, category, reason))

    connections = tuple(str(connection.name) for connection in getattr(structure, "connections", ()))
    return StructureAnalysis(path.resolve(), models, chains_by_model, tuple(components), protein_count, water_count,
                             tuple(sorted(altlocs)), tuple(sorted(insertion_codes)), tuple(zero_occupancy), connections,
                             tuple(protein_issues), tuple(protein_residues))


def atoms_for_residues(path: Path, model_index: int, keys: tuple[ResidueKey, ...]) -> tuple[AtomRecord, ...]:
    gemmi = _gemmi()
    structure = gemmi.read_structure(str(path))
    wanted = {(key.chain, key.number, key.insertion_code) for key in keys}
    atoms: list[AtomRecord] = []
    for chain in structure[model_index]:
        for residue in chain:
            identity = (chain.name, int(residue.seqid.num), str(residue.seqid.icode).strip())
            if identity not in wanted:
                continue
            atoms.extend(AtomRecord(atom.name.strip(), atom.element.name, float(atom.pos.x), float(atom.pos.y),
                                    float(atom.pos.z), float(atom.occ), str(atom.altloc).strip().replace("\x00", ""))
                         for atom in residue)
    if not atoms:
        raise ValueError("None of the selected binding-site residues were found")
    return tuple(atoms)
