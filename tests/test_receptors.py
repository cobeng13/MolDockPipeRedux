from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from moldockpipe.receptors.analysis import analysis_as_record, analyze_structure
from moldockpipe.receptors.box_calculation import center_from_atoms, envelope_box, radius_of_gyration_box
from moldockpipe.receptors.classification import classify_component
from moldockpipe.receptors.extraction import structure_as_pdb, write_cleaned_receptor, write_reference_ligand
from moldockpipe.receptors.models import AtomRecord, ComponentRole, ReceptorPreparationPlan, ResidueKey
from moldockpipe.project import ProjectRepository
from moldockpipe.receptors.preparation import _meeko_failure_message, prepare_receptor


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
    inventory = analysis_as_record(analysis, selected_model=0, included_chains=("A",))
    assert inventory["counts"]["protein_residues"] == 1
    assert inventory["protein_residues"] == ["ALA A:1"]
    assert {component["name"] for component in inventory["components"]} == {"LIG", "HOH", "ZN"}
    assert roles == {"LIG": ComponentRole.REFERENCE_LIGAND, "HOH": ComponentRole.REMOVE, "ZN": ComponentRole.RETAINED_ION}
    ligand = ResidueKey("A", 401, "", "LIG")
    plan = ReceptorPreparationPlan("test", "Test", source, 0, ("A",), ligand,
        (ResidueKey("A", 501, "", "HOH"),), (ResidueKey("A", 601, "", "ZN"),),
        (15, 16, 17), (20, 20, 20), "manual", center_method="reference_ligand_centroid",
        center_parameters={"reference_ligand": "LIG A:401"})
    assert plan.as_record()["center_method"] == "reference_ligand_centroid"
    pdb_text = structure_as_pdb(source)
    cleaned, reference = tmp_path / "cleaned.pdb", tmp_path / "reference.pdb"
    write_cleaned_receptor(pdb_text, plan, cleaned); write_reference_ligand(pdb_text, plan, reference)
    assert " LIG " not in cleaned.read_text(encoding="ascii")
    assert " HOH " not in cleaned.read_text(encoding="ascii")
    assert " ZN " in cleaned.read_text(encoding="ascii")
    assert reference.read_text(encoding="ascii").count("LIG") == 5


def test_meeko_failure_message_identifies_excess_bond_residues(tmp_path: Path) -> None:
    message = _meeko_failure_message(
        "matched with excess inter-residue bond(s): A:916\nmatched with excess inter-residue bond(s): A:942\n",
        tmp_path / "failure",
    )
    assert "A:916, A:942" in message
    assert "did not delete residues automatically" in message


def test_mm_cif_style_altloc_selection_uses_complete_conformer(tmp_path: Path) -> None:
    # This represents a residue whose only side-chain conformer is B.  The old
    # hard-coded A selection silently dropped CB before Meeko saw the PDB.
    source = tmp_path / "altloc.pdb"
    source.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N\n"
        "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00 20.00           C\n"
        "ATOM      3  C   ALA A   1       2.000   0.000   0.000  1.00 20.00           C\n"
        "ATOM      4  O   ALA A   1       3.000   0.000   0.000  1.00 20.00           O\n"
        "ATOM      5  CB BALA A   1       1.000   1.000   0.000  1.00 20.00           C\nEND\n",
        encoding="ascii",
    )
    analysis = analyze_structure(source)
    assert not analysis.protein_residue_issues
    plan = ReceptorPreparationPlan("test", "Test", source, 0, ("A",), None, (), (), (0, 0, 0), (20, 20, 20), "manual")
    cleaned = tmp_path / "cleaned.pdb"
    write_cleaned_receptor(structure_as_pdb(source), plan, cleaned)
    assert " CB " in cleaned.read_text(encoding="ascii")


def test_incomplete_standard_residue_is_reported(tmp_path: Path) -> None:
    source = tmp_path / "incomplete.pdb"
    source.write_text(PDB, encoding="ascii")
    issue = analyze_structure(source).protein_residue_issues[0]
    assert issue.key.label() == "ALA A:1"
    assert issue.missing_atoms == ("C", "CB", "O")


def test_receptor_preparation_persists_report_ready_provenance(tmp_path: Path, monkeypatch) -> None:
    repo = ProjectRepository.create(tmp_path / "project")
    source = tmp_path / "source.pdb"; source.write_text(PDB, encoding="ascii")
    plan = ReceptorPreparationPlan("rec-a", "Receptor A", source, 0, ("A",), None,
        (ResidueKey("A", 501, "", "HOH"), ResidueKey("A", 401, "", "LIG")),
        (ResidueKey("A", 601, "", "ZN"),), (15, 16, 17), (20, 21, 22),
        "manual", center_method="manual")

    def fake_meeko(command, *, cwd, **_kwargs):
        work = Path(cwd)
        (work / "receptor.pdbqt").write_text(
            "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00     0.000 N\n",
            encoding="utf-8")
        (work / "receptor.json").write_text("{}", encoding="utf-8")
        (work / "receptor.box.txt").write_text("box", encoding="utf-8")
        (work / "receptor_prepared.pdb").write_text(PDB, encoding="ascii")
        return SimpleNamespace(returncode=0, stdout="prepared", stderr="")

    monkeypatch.setattr("moldockpipe.receptors.preparation.subprocess.run", fake_meeko)
    folder = prepare_receptor(repo.root, plan)
    inventory = (folder / "structure_inventory.json").read_text(encoding="utf-8")
    report = (folder / "preparation_report.json").read_text(encoding="utf-8")
    with repo.connection() as conn:
        record = conn.execute("SELECT * FROM receptor_preparation_runs WHERE receptor_profile_id='rec-a'").fetchone()
    assert '"protein_residues": 1' in inventory
    assert '"method": "manual"' in report
    assert record["status"] == "completed"
    assert '"structure_inventory.json"' in record["artifacts_json"]
