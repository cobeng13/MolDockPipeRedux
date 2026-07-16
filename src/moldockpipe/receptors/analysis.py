from __future__ import annotations

from pathlib import Path

from .classification import WATER_NAMES, classify_component
from .models import AtomRecord, Component, ResidueKey, StructureAnalysis


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
                continue
            if water:
                water_count += 1
            polymer_like = amino
            category, role, reason = classify_component(name, atoms, polymer_like=polymer_like)
            components.append(Component(key, atoms, role, category, reason))

    connections = tuple(str(connection.name) for connection in getattr(structure, "connections", ()))
    return StructureAnalysis(path.resolve(), models, chains_by_model, tuple(components), protein_count, water_count,
                             tuple(sorted(altlocs)), tuple(sorted(insertion_codes)), tuple(zero_occupancy), connections)


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
