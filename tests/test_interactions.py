"""Geometry fixtures deliberately separate contacts from near-miss decoys."""
import math

import pytest

pytest.importorskip("rdkit")
from rdkit import Chem
from moldockpipe.services.interactions import analyze_interactions


def pdb_atom(serial, name, resn, xyz, element=None, chain="A", resi=10):
    element = element or name[0]
    return (f"ATOM  {serial:5d} {name:^4s} {resn:3s} {chain}{resi:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00 20.00          {element:>2s}\n")


def sdf(smiles, coordinates, explicit_h=False):
    mol = Chem.MolFromSmiles(smiles)
    if explicit_h:
        mol = Chem.AddHs(mol)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, coordinates[i])
    mol.AddConformer(conf)
    return Chem.MolToMolBlock(mol) + "\n$$$$\n"


def kinds(pdb, ligand):
    contacts, notes = analyze_interactions(pdb, ligand)
    assert "Geometric predictions" in notes
    return contacts


def hexagon(z=0, shift=0, tilt=0):
    return [(shift + 1.4 * math.cos(i * math.pi / 3),
             1.4 * math.sin(i * math.pi / 3) * math.cos(tilt),
             z + 1.4 * math.sin(i * math.pi / 3) * math.sin(tilt)) for i in range(6)]


def phenyl():
    return "".join(pdb_atom(i + 1, name, "PHE", xyz)
                   for i, (name, xyz) in enumerate(zip(("CG", "CD1", "CE1", "CZ", "CE2", "CD2"), hexagon())))


def test_hbond_candidate_needs_chemistry_distance_and_donor_geometry():
    protein = pdb_atom(1, "O", "ALA", (0, 0, 0)) + pdb_atom(2, "C", "ALA", (-1.2, 0, 0))
    ligand = sdf("CN", [(4.4, 0, 0), (3, 0, 0)])
    contacts = [c for c in kinds(protein, ligand) if c["kind"] == "hbond"]
    assert len(contacts) == 1
    assert "candidate" in contacts[0]["detail"]
    assert contacts[0]["distance"] == 3
    assert contacts[0]["residue"] == "ALA A:10"
    assert contacts[0]["start"] == {"x": 0, "y": 0, "z": 0}
    # Same distance, donor substituent points at the acceptor: reject.
    assert not [c for c in kinds(protein, sdf("CN", [(1.6, 0, 0), (3, 0, 0)])) if c["kind"] == "hbond"]
    assert not [c for c in kinds(protein, sdf("CN", [(5.4, 0, 0), (4, 0, 0)])) if c["kind"] == "hbond"]
    assert not [c for c in kinds(protein, sdf("CC", [(4.4, 0, 0), (3, 0, 0)])) if c["kind"] == "hbond"]


def test_explicit_hbond_uses_hydrogen_angle_and_distance():
    protein = pdb_atom(1, "O", "ALA", (0, 0, 0))
    # CN then three carbon H and two nitrogen H in RDKit atom order.
    coords = [(4.4, 0, 0), (3, 0, 0), (5.4, 0, 0), (4.4, 1, 0), (4.4, -1, 0), (2, 0, 0), (3, 1, 0)]
    contacts = [c for c in kinds(protein, sdf("CN", coords, True)) if c["kind"] == "hbond"]
    assert len(contacts) == 1 and "candidate" not in contacts[0]["detail"]
    coords[5] = (3, -1, 0)
    assert not [c for c in kinds(protein, sdf("CN", coords, True)) if c["kind"] == "hbond"]


@pytest.mark.parametrize("tilt,expected", [(0, "parallel"), (math.pi / 2, "T-shaped"), (math.pi / 4, None)])
def test_aromatic_stacking_checks_ring_normals(tilt, expected):
    contacts = [c for c in kinds(phenyl(), sdf("c1ccccc1", hexagon(z=4, tilt=tilt))) if c["kind"] == "pi_stack"]
    if expected:
        assert len(contacts) == 1 and expected in contacts[0]["detail"]
        assert contacts[0]["distance"] == 4
        assert len(contacts[0]["protein_serials"]) == 6
    else:
        assert not contacts


@pytest.mark.parametrize("smiles,coords", [
    ("C1CCCCC1", hexagon(z=3.5)),  # A planar hexagon alone is not aromatic.
    ("c1ccccc1", hexagon(z=3.5, shift=4)),
    ("c1ccccc1", hexagon(z=6)),
])
def test_stacking_rejects_nonaromatic_offset_and_distant_rings(smiles, coords):
    assert not [c for c in kinds(phenyl(), sdf(smiles, coords)) if c["kind"] == "pi_stack"]


def test_cation_pi_requires_actual_formal_charge_and_offset():
    contacts = kinds(phenyl(), sdf("[NH4+]", [(0, 0, 4)]))
    assert len([c for c in contacts if c["kind"] == "cation_pi"]) == 1
    assert not [c for c in kinds(phenyl(), sdf("N", [(0, 0, 4)])) if c["kind"] == "cation_pi"]
    assert not [c for c in kinds(phenyl(), sdf("[NH4+]", [(3, 0, 4)])) if c["kind"] == "cation_pi"]


def test_salt_bridge_both_charge_directions_and_acid_protonation():
    acid = pdb_atom(1, "OD1", "ASP", (-0.6, 0, 0)) + pdb_atom(2, "OD2", "ASP", (0.6, 0, 0))
    assert any(c["kind"] == "salt_bridge" for c in kinds(acid, sdf("[NH4+]", [(0, 0, 3.6)])))
    assert not any(c["kind"] == "salt_bridge" for c in kinds(acid, sdf("N", [(0, 0, 3.6)])))
    acid += pdb_atom(3, "HD1", "ASP", (-0.6, 0, 1), "H")
    assert not any(c["kind"] == "salt_bridge" for c in kinds(acid, sdf("[NH4+]", [(0, 0, 3.6)])))
    lysine = pdb_atom(1, "NZ", "LYS", (0, 0, 0))
    contacts = kinds(lysine, sdf("CC(=O)[O-]", [(0, 0, 4.5), (0, 0, 3.1), (1, 0, 2.5), (-1, 0, 2.5)]))
    assert any(c["kind"] == "salt_bridge" for c in contacts)
    assert any(c["kind"] == "cation_pi" for c in kinds(lysine, sdf("c1ccccc1", hexagon(z=4))))


def test_hydrophobic_contacts_deduplicate_per_residue_and_exclude_polar_carbon():
    protein = pdb_atom(1, "CD1", "LEU", (0, 0, 0)) + pdb_atom(2, "CG", "LEU", (-1.5, 0, 0))
    contacts = [c for c in kinds(protein, sdf("CCC", [(3, 0, 0), (3, 1.5, 0), (3, 3, 0)])) if c["kind"] == "hydrophobic"]
    assert len(contacts) == 1 and contacts[0]["distance"] == 3
    assert not [c for c in kinds(protein, sdf("C=O", [(3, 0, 0), (4.2, 0, 0)])) if c["kind"] == "hydrophobic"]


def test_unknown_residues_and_metals_are_not_assigned_contacts():
    protein = pdb_atom(1, "NZ", "UNK", (0, 0, 0)) + pdb_atom(2, "ZN", "ZN", (0, 0, 0), "Zn")
    assert kinds(protein, sdf("c1ccccc1", hexagon(z=4))) == []


def test_cancelled_and_invalid_chemistry_leave_coordinates_untouched():
    with pytest.raises(InterruptedError):
        analyze_interactions(phenyl(), sdf("N", [(0, 0, 3)]), cancelled=lambda: True)
    with pytest.raises(ValueError, match="chemistry"):
        analyze_interactions(phenyl(), "not an SDF")
