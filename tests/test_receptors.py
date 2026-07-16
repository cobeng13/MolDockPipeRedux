from __future__ import annotations

from pathlib import Path

import pytest

from moldockpipe.receptors.analysis import analyze_structure
from moldockpipe.receptors.box_calculation import center_from_atoms, envelope_box, radius_of_gyration_box
from moldockpipe.receptors.classification import classify_component
from moldockpipe.receptors.extraction import structure_as_pdb, write_cleaned_receptor, write_reference_ligand
from moldockpipe.receptors.models import AtomRecord, ComponentRole, ReceptorPreparationPlan, ResidueKey


PDB = """\
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00 20.00           C
HETATM    3  C1  LIG A 401      14.000  15.000  16.000  1.00 20.00           C
HETATM    4  C2  LIG A 401      16.000  17.000  18.000  1.00 20.00           C
HETATM    5  O1  LIG A 401      15.000  16.000  17.000  1.00 20.00           O
HETATM    6  N1  LIG A 401      15.000  15.000  17.000  1.00 20.00           N
HETATM    7  C3  LIG A 401      16.000  15.000  17.000  1.00 20.00           C
HETATM    8  O   HOH A 501      20.000  20.000  20.000  1.00 20.00           O
HETATM    9 ZN    ZN A 601      12.000  12.000  12.000  1.00 20.00          ZN
END
"""


def atom(x: float, y: float, z: float, element: str = "C") -> AtomRecord:
    return AtomRecord("C", element, x, y, z, 1.0)


def test_box_calculations_use_heavy_atoms() -> None:
    atoms = (atom(0, 0, 0), atom(2, 4, 6), atom(99, 99, 99, "H"))
    assert center_from_atoms(atoms) == (1, 2, 3)
    center, size = envelope_box(atoms, 8)
    assert center == (1, 2, 3)
    assert size == (18, 20, 22)
    _, cubic = radius_of_gyration_box(atoms)
    assert cubic[0] == pytest.approx(cubic[1]) == pytest.approx(cubic[2])


def test_component_classification_is_conservative() -> None:
    assert classify_component("HOH", (atom(0, 0, 0, "O"),))[1] == ComponentRole.REMOVE
    assert classify_component("ZN", (atom(0, 0, 0, "Zn"),))[1] == ComponentRole.RETAINED_ION
    assert classify_component("HEM", tuple(atom(i, 0, 0) for i in range(5)))[1] == ComponentRole.RETAINED_COFACTOR
    assert classify_component("UNK", (atom(0, 0, 0),))[1] == ComponentRole.UNRESOLVED


def test_analysis_and_extraction(tmp_path: Path) -> None:
    source = tmp_path / "source.pdb"; source.write_text(PDB, encoding="ascii")
    analysis = analyze_structure(source)
    roles = {component.key.name: component.suggested_role for component in analysis.components}
    assert analysis.protein_residue_count == 1
    assert roles == {"LIG": ComponentRole.REFERENCE_LIGAND, "HOH": ComponentRole.REMOVE, "ZN": ComponentRole.RETAINED_ION}
    ligand = ResidueKey("A", 401, "", "LIG")
    plan = ReceptorPreparationPlan("test", "Test", source, 0, ("A",), ligand,
        (ResidueKey("A", 501, "", "HOH"),), (ResidueKey("A", 601, "", "ZN"),),
        (15, 16, 17), (20, 20, 20), "manual")
    pdb_text = structure_as_pdb(source)
    cleaned, reference = tmp_path / "cleaned.pdb", tmp_path / "reference.pdb"
    write_cleaned_receptor(pdb_text, plan, cleaned); write_reference_ligand(pdb_text, plan, reference)
    assert " LIG " not in cleaned.read_text(encoding="ascii")
    assert " HOH " not in cleaned.read_text(encoding="ascii")
    assert " ZN " in cleaned.read_text(encoding="ascii")
    assert reference.read_text(encoding="ascii").count("LIG") == 5
