from pathlib import Path

from moldockpipe.services.postdock import PostDockService, receptor_export_filename


PDB_ATOM = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nEND\n"
PDBQT_ATOM = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00    -0.123 C\n"
MOL2 = """@<TRIPOS>MOLECULE
reference_ligand
 2 1 0 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.3       1 LIG       0.0000
      2 O1          1.2000    0.0000    0.0000 O.3       1 LIG       0.0000
@<TRIPOS>BOND
     1    1    2 1
"""


def test_postdock_exports_meeko_prepared_pdb_not_cleaned_input(tmp_path: Path) -> None:
    cleaned = tmp_path / "receptor_cleaned.pdb"
    prepared = tmp_path / "receptor_prepared.pdb"
    pdbqt = tmp_path / "receptor.pdbqt"
    output = tmp_path / "bundle" / "COX-2_4M11_prepared.pdb"
    cleaned.write_text(PDB_ATOM.replace("1.000", "9.000"), encoding="ascii")
    prepared.write_text(PDB_ATOM, encoding="ascii")
    pdbqt.write_text(PDBQT_ATOM, encoding="ascii")

    source = PostDockService().export_receptor_pdb(pdbqt, output, prepared)

    assert source == "meeko_prepared_pdb"
    assert output.read_text(encoding="ascii") == PDB_ATOM
    assert "9.000" not in output.read_text(encoding="ascii")
    assert receptor_export_filename("COX-2 4M11") == "COX-2_4M11_prepared.pdb"


def test_postdock_can_derive_pdb_from_imported_prepared_pdbqt(tmp_path: Path) -> None:
    pdbqt = tmp_path / "receptor.pdbqt"
    output = tmp_path / "bundle" / "Imported_prepared.pdb"
    pdbqt.write_text(PDBQT_ATOM, encoding="ascii")

    source = PostDockService().export_receptor_pdb(pdbqt, output)

    text = output.read_text(encoding="ascii")
    assert source == "prepared_pdbqt_conversion"
    assert "Generated from the prepared receptor PDBQT" in text
    assert "ATOM" in text and text.rstrip().endswith("END")


def test_postdock_copies_exact_validation_mol2_standard(tmp_path: Path) -> None:
    source = tmp_path / "redocking" / "dockrmsd" / "reference_ligand_heavy.mol2"
    output = tmp_path / "bundle" / "reference_ligand_heavy.mol2"
    source.parent.mkdir(parents=True)
    source.write_text(MOL2, encoding="utf-8")

    PostDockService().export_validation_standard(source, output)

    assert output.read_bytes() == source.read_bytes()
