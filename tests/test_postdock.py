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


def test_split_stages_long_paths_and_preserves_names(tmp_path, monkeypatch):
    from types import SimpleNamespace

    source = tmp_path / ('input_' + 'x' * 80) / 'raw.pdbqt'
    source.parent.mkdir()
    source.write_text('MODEL 1\nENDMDL\n', encoding='ascii')
    output = tmp_path / ('a' * 90) / ('b' * 90) / ('c' * 50)
    prefix = 'Quercetin7Rhamnoside_Quercetin7Rhamnoside-0eab8f797a06dd29_pose_'
    staging = []

    def run(args, **kwargs):
        work = Path(kwargs['cwd'])
        staging.append(work)
        assert args[2:] == ['input.pdbqt', '--ligand', 'pose_']
        assert (work / 'input.pdbqt').read_bytes() == source.read_bytes()
        for number in (1, 2):
            (work / f'pose_{number}.pdbqt').write_text(f'pose {number}', encoding='ascii')
        return SimpleNamespace(returncode=0, stdout='', stderr='')

    monkeypatch.setattr('moldockpipe.services.postdock.subprocess.run', run)
    outputs = PostDockService().split(Path('vina_split.exe'), source, output, prefix)
    assert outputs == [output / f'{prefix}{number}.pdbqt' for number in (1, 2)]
    assert [p.read_text() for p in outputs] == ['pose 1', 'pose 2']
    assert not staging[0].exists()


def test_sdf_stages_long_paths_and_keeps_log(tmp_path, monkeypatch):
    from types import SimpleNamespace

    source = tmp_path / 'pose.pdbqt'
    source.write_text('pose', encoding='ascii')
    output = tmp_path / ('a' * 90) / ('b' * 90) / ('c' * 50) / 'pose.sdf'
    log = output.with_suffix('.log.txt')
    staging = []

    def run(args, **kwargs):
        work = Path(kwargs['cwd'])
        staging.append(work)
        assert args[-3:] == ['input.pdbqt', '-s', 'output.sdf']
        assert (work / 'input.pdbqt').read_bytes() == source.read_bytes()
        (work / 'output.sdf').write_text('molecule\n$$$$\n', encoding='ascii')
        return SimpleNamespace(returncode=0, stdout='exported', stderr='warning')

    monkeypatch.setattr('moldockpipe.services.postdock.subprocess.run', run)
    PostDockService().export_sdf(source, output, log)
    assert output.read_text() == 'molecule\n$$$$\n'
    assert log.read_text() == 'exported\nwarning'
    assert not staging[0].exists()


def test_split_failure_cleans_staging_without_publishing(tmp_path, monkeypatch):
    import pytest
    from types import SimpleNamespace

    source = tmp_path / 'raw.pdbqt'
    source.write_text('pose', encoding='ascii')
    output = tmp_path / 'exports'
    staging = []

    def run(args, **kwargs):
        work = Path(kwargs['cwd'])
        staging.append(work)
        (work / 'pose_1.pdbqt').write_text('incomplete', encoding='ascii')
        return SimpleNamespace(returncode=1, stdout='', stderr='split failed')

    monkeypatch.setattr('moldockpipe.services.postdock.subprocess.run', run)
    with pytest.raises(RuntimeError, match='split failed'):
        PostDockService().split(Path('vina_split.exe'), source, output, 'pose_')
    assert not output.exists()
    assert not staging[0].exists()
