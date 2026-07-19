from __future__ import annotations

import json
from collections import Counter
from pathlib import Path


def _chem():
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise RuntimeError("RDKit is required for reference-ligand chemistry reconstruction") from exc
    return Chem


def read_sdf(path: Path, *, remove_hydrogens: bool = False):
    Chem = _chem()
    supplier = Chem.SDMolSupplier(str(path), removeHs=remove_hydrogens, sanitize=True)
    molecule = next((mol for mol in supplier if mol is not None), None)
    if molecule is None or molecule.GetNumConformers() == 0:
        raise ValueError(f"SDF is unreadable or has no 3D molecule: {path}")
    return molecule


def molecule_from_crystal_pdb(path: Path, template_sdf: Path | None = None):
    Chem = _chem()
    molecule = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=template_sdf is None, proximityBonding=True)
    if molecule is None or molecule.GetNumAtoms() == 0 or molecule.GetNumConformers() == 0:
        raise ValueError("RDKit could not infer complete ligand chemistry from the crystallographic PDB")
    if template_sdf is not None:
        from rdkit.Chem import AllChem
        template = read_sdf(template_sdf, remove_hydrogens=True); crystal = Chem.RemoveHs(molecule, sanitize=False)
        expected = Counter(atom.GetSymbol() for atom in template.GetAtoms()); observed = Counter(atom.GetSymbol() for atom in crystal.GetAtoms())
        if template.GetNumAtoms() != crystal.GetNumAtoms() or expected != observed:
            raise ValueError("Chemical template does not map completely to crystallographic atoms. "
                             f"Template atoms: {template.GetNumAtoms()}; crystallographic atoms: {crystal.GetNumAtoms()}")
        try: molecule = AllChem.AssignBondOrdersFromTemplate(template, crystal)
        except Exception as exc: raise ValueError("Chemical template connectivity could not be mapped completely to crystallographic coordinates") from exc
        Chem.SanitizeMol(molecule); Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    if molecule.GetNumBonds() == 0 and molecule.GetNumHeavyAtoms() > 1:
        raise ValueError("Reference ligand has multiple atoms but no inferred bonds")
    return molecule


def atom_name(atom, index: int) -> str:
    info = atom.GetPDBResidueInfo()
    if info and info.GetName().strip():
        return info.GetName().strip()
    if atom.HasProp("_TriposAtomName"):
        return atom.GetProp("_TriposAtomName")
    return f"{atom.GetSymbol()}{index + 1}"


def chemistry_summary(molecule) -> dict[str, object]:
    Chem = _chem(); heavy = Chem.RemoveHs(molecule, sanitize=True)
    return {
        "heavy_atom_count": heavy.GetNumAtoms(),
        "elements": dict(sorted(Counter(atom.GetSymbol() for atom in heavy.GetAtoms()).items())),
        "formal_charge": sum(atom.GetFormalCharge() for atom in molecule.GetAtoms()),
        "canonical_smiles": Chem.MolToSmiles(heavy, isomericSmiles=True),
    }


def write_sdf(molecule, destination: Path) -> None:
    Chem = _chem(); destination.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(destination)); writer.write(molecule); writer.close()
    if not destination.is_file() or destination.stat().st_size == 0:
        raise RuntimeError(f"Failed to write SDF: {destination}")


def _mol2_atom_type(atom) -> str:
    from rdkit.Chem.rdchem import HybridizationType
    symbol = atom.GetSymbol()
    if atom.GetIsAromatic(): return f"{symbol}.ar"
    if symbol == "C": return "C.2" if atom.GetHybridization() in {HybridizationType.SP, HybridizationType.SP2} else "C.3"
    if symbol == "N": return "N.2" if atom.GetHybridization() == HybridizationType.SP2 else "N.3"
    if symbol == "O": return "O.2" if atom.GetHybridization() == HybridizationType.SP2 else "O.3"
    return symbol


def write_mol2(molecule, destination: Path, *, heavy_only: bool = True, name: str = "reference_ligand") -> None:
    Chem = _chem()
    mol = Chem.RemoveHs(molecule, sanitize=True) if heavy_only else Chem.Mol(molecule)
    if mol.GetNumConformers() == 0: raise ValueError("MOL2 export requires 3D coordinates")
    conf = mol.GetConformer(); destination.parent.mkdir(parents=True, exist_ok=True)
    lines = ["@<TRIPOS>MOLECULE", name, f"{mol.GetNumAtoms()} {mol.GetNumBonds()} 1 0 0", "SMALL", "USER_CHARGES", "",
             "@<TRIPOS>ATOM"]
    for index, atom in enumerate(mol.GetAtoms(), 1):
        point = conf.GetAtomPosition(index - 1)
        lines.append(f"{index:7d} {atom_name(atom, index-1):<8} {point.x:10.4f} {point.y:10.4f} {point.z:10.4f} "
                     f"{_mol2_atom_type(atom):<6} 1 LIG {float(atom.GetFormalCharge()):9.4f}")
    lines.append("@<TRIPOS>BOND")
    for index, bond in enumerate(mol.GetBonds(), 1):
        kind = "ar" if bond.GetIsAromatic() else {1.0: "1", 2.0: "2", 3.0: "3"}.get(bond.GetBondTypeAsDouble(), "1")
        lines.append(f"{index:6d} {bond.GetBeginAtomIdx()+1:5d} {bond.GetEndAtomIdx()+1:5d} {kind}")
    lines.extend(("@<TRIPOS>SUBSTRUCTURE", "     1 LIG         1 GROUP", ""))
    destination.write_text("\n".join(lines), encoding="utf-8")


def validate_same_heavy_graph(reference, pose) -> None:
    match_heavy_atom_graph(reference, pose)


def match_heavy_atom_graph(reference, pose) -> tuple[int, ...]:
    """Return a pose-to-reference atom mapping for the same heavy-atom graph.

    Meeko's PDBQT round trip may reorder atoms and can assign stereochemical
    tags that were absent from crystallographic PDB/mmCIF coordinates.  Neither
    changes molecular identity for the external RMSD comparison, so graph
    matching deliberately ignores atom order and stereochemical labels while
    retaining element and bond-order/aromaticity checks.
    """
    Chem = _chem()
    ref_mol = Chem.RemoveHs(reference, sanitize=True)
    pose_mol = Chem.RemoveHs(pose, sanitize=True)
    ref, other = chemistry_summary(ref_mol), chemistry_summary(pose_mol)
    if ref["heavy_atom_count"] != other["heavy_atom_count"] or ref["elements"] != other["elements"]:
        raise ValueError("The redocked pose molecular graph does not match the reference ligand. "
                         f"Reference heavy atoms: {ref['heavy_atom_count']}; pose heavy atoms: {other['heavy_atom_count']}")
    # The match tuple is indexed by reference atom order and contains the
    # matching atom index in the pose. useChirality=False is essential here:
    # reference structures normally do not encode stereochemistry.
    mapping = pose_mol.GetSubstructMatch(ref_mol, useChirality=False)
    if len(mapping) != ref_mol.GetNumAtoms():
        raise ValueError("The redocked pose bond graph does not match the reference ligand")
    for reference_index, pose_index in enumerate(mapping):
        if ref_mol.GetAtomWithIdx(reference_index).GetFormalCharge() != pose_mol.GetAtomWithIdx(pose_index).GetFormalCharge():
            raise ValueError("The redocked pose formal-charge graph does not match the reference ligand")
    return tuple(mapping)


def pose_in_reference_order(reference, pose):
    """Return a heavy-atom pose ordered and named exactly like its reference."""
    Chem = _chem()
    reference_heavy = Chem.RemoveHs(reference, sanitize=True)
    pose_heavy = Chem.RemoveHs(pose, sanitize=True)
    mapping = match_heavy_atom_graph(reference_heavy, pose_heavy)
    aligned = Chem.RenumberAtoms(pose_heavy, list(mapping))
    for index, atom in enumerate(aligned.GetAtoms()):
        # DockRMSD and manual inspection can now use the same stable labels in
        # every generated MOL2, even when Meeko changed the PDBQT atom order.
        atom.SetProp("_TriposAtomName", atom_name(reference_heavy.GetAtomWithIdx(index), index))
        atom.SetIntProp("MOLDOCK_REFERENCE_ATOM_INDEX", index + 1)
    return aligned


def create_reference_bundle(pdb_path: Path, folder: Path, *, template_sdf: Path | None = None,
                            chemistry_source: str | None = None) -> dict[str, object]:
    """Create immutable coordinate-preserving chemistry files; raises on incomplete inference."""
    molecule = molecule_from_crystal_pdb(pdb_path, template_sdf)
    folder.mkdir(parents=True, exist_ok=True)
    bundled_pdb = folder / "reference_ligand.pdb"
    if pdb_path.resolve() != bundled_pdb.resolve(): bundled_pdb.write_bytes(pdb_path.read_bytes())
    sdf = folder / "reference_ligand.sdf"; mol2 = folder / "reference_ligand.mol2"
    write_sdf(molecule, sdf); write_mol2(molecule, mol2, heavy_only=True)
    mapping = {
        "chemistry_source": chemistry_source or ("user_supplied_sdf" if template_sdf else "rdkit_pdb_inference"),
        "atoms": [{"index": i, "name": atom_name(atom, i), "element": atom.GetSymbol()} for i, atom in enumerate(molecule.GetAtoms())],
        **chemistry_summary(molecule),
    }
    (folder / "chemistry_mapping.json").write_text(json.dumps(mapping, indent=2), encoding="utf-8")
    return mapping


def hydrogenated_copy(reference_sdf: Path, destination: Path) -> None:
    Chem = _chem(); molecule = read_sdf(reference_sdf)
    before = molecule.GetConformer(); coordinates = [(before.GetAtomPosition(i).x, before.GetAtomPosition(i).y, before.GetAtomPosition(i).z)
                                                      for i, atom in enumerate(molecule.GetAtoms()) if atom.GetAtomicNum() > 1]
    hydrogenated = Chem.AddHs(molecule, addCoords=True)
    after = hydrogenated.GetConformer(); actual = [(after.GetAtomPosition(i).x, after.GetAtomPosition(i).y, after.GetAtomPosition(i).z)
                                                   for i, atom in enumerate(hydrogenated.GetAtoms()) if atom.GetAtomicNum() > 1]
    if coordinates != actual: raise RuntimeError("Hydrogenation changed crystallographic heavy-atom coordinates")
    write_sdf(hydrogenated, destination)
