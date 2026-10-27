import csv
import json
from pathlib import Path

import pytest

from moldockpipe.project import ProjectRepository
from moldockpipe.fingerprints import file_sha256
from moldockpipe.services.interaction_analysis import (
    InteractionAnalysisService, ProlifEngine, compare_fingerprints, contact_pair,
    heavy_coordinates, pose_mapping, same_coordinates, sdf_records, validate_molecule,
)

Chem = pytest.importorskip("rdkit.Chem")


def molecule(smiles="CN", points=None):
    mol = Chem.MolFromSmiles(smiles)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, (points or [(i * 1.4, 0, 0) for i in range(mol.GetNumAtoms())])[i])
    mol.AddConformer(conf)
    return mol


def write_sdf(path, molecules):
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(path))
    for mol in molecules:
        writer.write(mol)
    writer.close()


def pdb_line(i, name, xyz, *, element="C", residue="ALA", number=10, chain="A", insertion=""):
    return (f"ATOM  {i:5d} {name:^4s} {residue:3s} {chain}{number:4d}{insertion or ' '}   "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00 20.00          {element:>2s}\n")


def raw_pose(mol):
    return "".join(pdb_line(i+1, symbol, xyz, element=symbol)[:76] + f"{symbol}\n"
                   for i, (symbol, xyz) in enumerate(heavy_coordinates(mol)))


@pytest.fixture
def project(tmp_path):
    repo = ProjectRepository.create(tmp_path / "project")
    root = repo.root
    receptor_dir = root / "inputs/receptors/r"
    receptor_dir.mkdir(parents=True)
    ligand = molecule()
    clean = pdb_line(1, "CA", (10, 0, 0))
    crystal = "".join(pdb_line(i+2, symbol+str(i), xyz, element=symbol, residue="LIG", number=20)
                      for i, (symbol, xyz) in enumerate(heavy_coordinates(ligand)))
    (receptor_dir / "source.pdb").write_text(clean + crystal + "END\n")
    for name in ("receptor_cleaned.pdb", "receptor_prepared.pdb", "receptor.pdbqt"):
        (receptor_dir / name).write_text(clean + "END\n")
    (receptor_dir / "reference.pdb").write_text(crystal + "END\n")
    write_sdf(receptor_dir / "reference.sdf", [ligand])
    report = {key: file_sha256(receptor_dir / name) for key, name in (
        ("source_sha256", "source.pdb"), ("cleaned_receptor_sha256", "receptor_cleaned.pdb"),
        ("prepared_receptor_sha256", "receptor.pdbqt"))}
    (receptor_dir / "preparation_report.json").write_text(json.dumps(report))
    profile = {"id": "r", "name": "Receptor", "enabled": True, "receptor": "inputs/receptors/r/receptor.pdbqt",
               "reference_ligand": {"identity": "LIG A20", "pdb": "inputs/receptors/r/reference.pdb", "sdf": "inputs/receptors/r/reference.sdf"}}
    repo.save_receptor_profiles([profile])
    raw = root / "artifacts/docking/raw.pdbqt"
    raw.write_text("MODEL 1\n" + raw_pose(ligand) + "ENDMDL\nMODEL 2\n" + raw_pose(molecule(points=[(0, 2, 0), (1.4, 2, 0)])) + "ENDMDL\n")
    artifact = repo.add_artifact(raw, "vina_output_pdbqt", "vina")
    with repo.connection() as conn:
        conn.execute("INSERT INTO parent_ligands(parent_id,source_smiles,parse_status,created_at) VALUES ('parent','CN','valid','now')")
        conn.execute("""INSERT INTO molecular_states(state_id,parent_id,state_smiles,state_isomeric_smiles,state_structure_hash,generation_fingerprint,status,created_at,formal_charge)
                     VALUES ('state','parent','CN','CN','hash','fp','ready','now',0)""")
        conn.execute("""INSERT INTO docking_runs(run_id,state_id,receptor_hash,settings_fingerprint,status,raw_output_artifact_id,is_current,receptor_profile_id,command_json)
                     VALUES ('run','state',?,'fp','completed',?,1,'r','{}')""", (report["prepared_receptor_sha256"], artifact))
        conn.executemany("INSERT INTO docking_poses(pose_id,run_id,mode_index,affinity) VALUES (?,'run',?,?)",
                         [("pose1", 1, -7), ("pose2", 2, -6)])
    output = root / "For_PostDocking/r/SDF/parent/state/run/all.sdf"
    write_sdf(output, [ligand, molecule(points=[(0, 2, 0), (1.4, 2, 0)])])
    return repo, output


CONTACT = dict(protein_chain="A", protein_residue="ALA", residue_number=10,
               insertion_code="", interaction_type="HBDonor", distance_angstrom=3.0)


class FakeEngine:
    calls = 0
    loads = 0
    residues = {1: ("A", "ALA", 10, "")}
    excluded = []

    def __init__(self, path):
        FakeEngine.loads += 1

    def analyze(self, mol):
        FakeEngine.calls += 1
        return [CONTACT.copy(), CONTACT.copy()]


def service(repo, **kwargs):
    return InteractionAnalysisService(repo, engine_factory=FakeEngine, **kwargs)


def reports(repo):
    with repo.connection() as conn:
        row = conn.execute("SELECT * FROM interaction_runs WHERE output_path IS NOT NULL ORDER BY finished_at DESC, rowid DESC LIMIT 1").fetchone()
    folder = repo.root / row["output_path"]
    return folder, list(csv.DictReader((folder / "interaction_summary.csv").open(encoding="utf-8-sig")))


def test_exact_pairs_and_residue_conservation():
    a = ("A", "ASP", 10, "", "HBDonor")
    b = ("B", "ASP", 10, "", "HBDonor")
    c = ("A", "ASP", 10, "A", "HBDonor")
    different_type = (*a[:4], "Hydrophobic")
    metrics = compare_fingerprints({a, b, different_type}, {a, c})
    assert metrics["conserved_interaction_pairs"] == 1
    assert metrics["lost_reference_pairs"] == 1
    assert metrics["new_interaction_pairs"] == 2
    assert metrics["jaccard_similarity"] == .25
    assert metrics["reference_retention_percent"] == 50
    assert compare_fingerprints({different_type}, {a})["residue_retention_percent"] == 100
    assert compare_fingerprints(set(), set())["jaccard_similarity"] is None
    assert compare_fingerprints({a}, set())["reference_retention_percent"] is None
    assert compare_fingerprints(set(), {a})["reference_retention_percent"] == 0


def test_streaming_invalid_empty_and_no_conformer(tmp_path):
    path = tmp_path / "multi.sdf"
    write_sdf(path, [molecule(), molecule()])
    with path.open("a") as handle:
        handle.write("bad molecule\n$$$$\n")
    rows = list(sdf_records(path))
    assert len(rows) == 3 and rows[2][1] is None
    with pytest.raises(ValueError, match="conformer"):
        validate_molecule(Chem.MolFromSmiles("CN"))
    with pytest.raises(FileNotFoundError):
        list(sdf_records(tmp_path / "missing.sdf"))
    path.write_text("")
    assert list(sdf_records(path)) == [(1, None)]
    mol = molecule()
    mol.GetConformer().SetAtomPosition(0, (float("nan"), 0, 0))
    with pytest.raises(ValueError, match="Non-finite"):
        validate_molecule(mol)


def test_coordinates_order_and_rounding():
    atoms = [("C", (0, 0, 0)), ("N", (1.4004, 0, 0))]
    assert same_coordinates(atoms, [("N", (1.4, 0, 0)), ("C", (0, 0, 0))])
    assert not same_coordinates(atoms, [("N", (1.41, 0, 0)), ("C", (0, 0, 0))])


def test_batch_csv_reuse_and_force(project):
    repo, path = project
    FakeEngine.calls = FakeEngine.loads = 0
    before = file_sha256(path)
    assert service(repo).run(profile_id="r") == (2, 0)
    assert FakeEngine.calls == 3 and FakeEngine.loads == 1
    folder, rows = reports(repo)
    assert [r["pose_id"] for r in rows] == ["pose1", "pose2"]
    assert [float(r["affinity"]) for r in rows] == [-7, -6]
    assert all(r["total_interaction_pairs"] == "1" for r in rows)
    detail = list(csv.DictReader((folder / "interactions_detailed.csv").open(encoding="utf-8-sig")))
    assert len(detail) == 6  # detailed duplicates retained, fingerprint deduplicated
    assert service(repo).run(profile_id="r") == (2, 0)
    assert FakeEngine.calls == 3
    assert service(repo).run(profile_id="r", force=True) == (2, 0)
    assert reports(repo)[0] != folder
    assert file_sha256(path) == before


def test_failed_record_is_not_zero_and_changed_inputs_invalidate(project):
    repo, path = project
    service(repo).run(profile_id="r")
    with path.open("a") as handle:
        handle.write("bad\n$$$$\n")
    assert service(repo).run(profile_id="r") == (2, 1)
    folder, rows = reports(repo)
    assert rows[2]["analysis_status"] == "failed"
    assert rows[2]["total_interaction_pairs"] == "NA"
    matrix = list(csv.DictReader((folder / "interaction_comparison.csv").open(encoding="utf-8-sig")))
    assert matrix[-1]["candidate_present"] == "NA" and matrix[-1]["comparison"] == "unknown"


@pytest.mark.parametrize("reference_empty", [True, False])
def test_valid_zero_contacts_and_empty_reference(project, reference_empty):
    repo, _ = project
    class EmptyEngine(FakeEngine):
        def analyze(self, mol):
            return [] if reference_empty or mol.GetConformer().GetAtomPosition(0).y else [CONTACT]
    analyzer = InteractionAnalysisService(repo, engine_factory=EmptyEngine)
    assert analyzer.run(profile_id="r") == (2, 0)
    _, rows = reports(repo)
    assert rows[1]["analysis_status"] == "success"
    assert rows[1]["total_interaction_pairs"] == "0"
    assert rows[1]["reference_retention_percent"] == ("NA" if reference_empty else "0.0")


def test_multiple_files_and_ambiguous_pose_mapping(project):
    repo, path = project
    write_sdf(path.with_name("extra.sdf"), [molecule()])
    with repo.connection() as conn:
        conn.execute("UPDATE docking_poses SET mode_index=3 WHERE pose_id='pose2'")
    assert service(repo).run(profile_id="r") == (3, 0)
    _, rows = reports(repo)
    ambiguous = next(r for r in rows if "ambiguous" in r["mapping_status"])
    assert ambiguous["affinity"] == "NA"


@pytest.mark.parametrize("file", ["reference.sdf", "receptor_prepared.pdb", "preparation_report.json"])
def test_missing_required_inputs_preserve_reports(project, file):
    repo, _ = project
    service(repo).run(profile_id="r")
    folder, _ = reports(repo)
    pointer = (folder.parent / "latest.json").read_bytes()
    (repo.root / "inputs/receptors/r" / file).unlink()
    with pytest.raises(FileNotFoundError):
        service(repo).run(profile_id="r")
    assert (folder.parent / "latest.json").read_bytes() == pointer


def test_coordinate_mismatch_and_cancellation_preserve_output(project):
    repo, path = project
    service(repo).run(profile_id="r")
    folder, _ = reports(repo)
    pointer = (folder.parent / "latest.json").read_bytes()
    write_sdf(path, [molecule(points=[(50, 0, 0), (51.4, 0, 0)])])
    assert service(repo).run(profile_id="r") == (0, 1)
    assert (folder.parent / "latest.json").read_bytes() == pointer
    with pytest.raises(InterruptedError):
        service(repo, cancelled=lambda: True).run(profile_id="r")


def test_fatal_calculation_preserves_reports(project):
    repo, _ = project
    service(repo).run(profile_id="r")
    folder, _ = reports(repo)
    class BrokenEngine(FakeEngine):
        def analyze(self, mol):
            raise RuntimeError("reference failure")
    with pytest.raises(RuntimeError, match="reference failure"):
        InteractionAnalysisService(repo, engine_factory=BrokenEngine).run(profile_id="r", force=True)
    assert folder.is_dir()
    assert not list(folder.parent.glob(".pending-*"))


def test_real_prolif_hbond_geometry_directionality_and_metals(tmp_path):
    # Full chemically valid alanine, with explicit hydrogens and PDB atom names.
    from rdkit.Chem import AllChem
    peptide = Chem.AddHs(Chem.MolFromSequence("A"), addResidueInfo=True)
    assert AllChem.EmbedMolecule(peptide, randomSeed=42) == 0
    path = tmp_path / "protein.pdb"
    Chem.MolToPDBFile(peptide, str(path))
    engine = ProlifEngine(path)
    # Select the carboxyl carbonyl O; put a donor N at 2.9 A with a linear N-H-O.
    oxygen = next(a for a in peptide.GetAtoms() if a.GetSymbol() == "O" and
                  any(b.GetBondTypeAsDouble() == 2 for b in a.GetBonds()))
    point = peptide.GetConformer().GetAtomPosition(oxygen.GetIdx())
    donor = Chem.AddHs(Chem.MolFromSmiles("N"))
    conf = Chem.Conformer(donor.GetNumAtoms())
    for i, delta in enumerate([(2.9, 0, 0), (1.9, 0, 0), (3.3, 1, 0), (3.3, -1, 0)]):
        conf.SetAtomPosition(i, tuple(point[j]+delta[j] for j in range(3)))
    donor.AddConformer(conf)
    contacts = engine.analyze(donor)
    assert any(c["interaction_type"] == "HBDonor" and c["angle_degrees"] > 170 for c in contacts)
    import numpy as np
    nitrogen = next(a for a in peptide.GetAtoms() if a.GetSymbol() == "N")
    hydrogen = next(a for a in nitrogen.GetNeighbors() if a.GetSymbol() == "H")
    n = np.array(peptide.GetConformer().GetAtomPosition(nitrogen.GetIdx()))
    h = np.array(peptide.GetConformer().GetAtomPosition(hydrogen.GetIdx()))
    axis = (h-n) / np.linalg.norm(h-n)
    acceptor_point = n + 2.9 * axis
    acceptor = molecule("C=O", [tuple(acceptor_point+(0, 1.2, 0)), tuple(acceptor_point)])
    assert any(c["interaction_type"] == "HBAcceptor" for c in engine.analyze(acceptor))
    # Move all atoms away: a valid zero-contact ligand, not a calculation failure.
    for i in range(donor.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        donor.GetConformer().SetAtomPosition(i, (p.x+100, p.y, p.z))
    assert engine.analyze(donor) == []
    assert set(engine.fp.interactions) >= {"MetalAcceptor", "XBDonor", "HBAcceptor"}


def test_real_metal_halogen_and_insertion_codes(tmp_path):
    import numpy as np
    from rdkit.Chem import AllChem
    peptide = Chem.AddHs(Chem.MolFromSequence("A"), addResidueInfo=True)
    AllChem.EmbedMolecule(peptide, randomSeed=42)
    for atom in peptide.GetAtoms():
        atom.GetPDBResidueInfo().SetResidueNumber(10)
        atom.GetPDBResidueInfo().SetInsertionCode("A")
    path = tmp_path / "protein.pdb"
    text = Chem.MolToPDBBlock(peptide).split("CONECT")[0].split("END")[0]
    text += pdb_line(999, "ZN", (30, 0, 0), element="Zn", residue="ZN", number=11).replace("ATOM  ", "HETATM")
    # Same residue number in another chain is separately represented.
    other = Chem.Mol(peptide)
    for atom in other.GetAtoms():
        atom.GetPDBResidueInfo().SetChainId("B")
        atom.GetPDBResidueInfo().SetInsertionCode("")
        p = other.GetConformer().GetAtomPosition(atom.GetIdx())
        other.GetConformer().SetAtomPosition(atom.GetIdx(), (p.x, p.y+30, p.z))
    text += Chem.MolToPDBBlock(other).split("CONECT")[0].split("END")[0] + "END\n"
    path.write_text(text)
    engine = ProlifEngine(path)
    assert len(engine.protein.residues) == 3
    assert ("A", "ALA", 10, "A") in engine.residues.values()
    contacts = engine.analyze(molecule("N", [(32, 0, 0)]))
    assert any(c["interaction_type"] == "MetalAcceptor" and c["protein_residue"] == "ZN" for c in contacts)
    assert not any(c["interaction_type"] == "MetalAcceptor" for c in engine.analyze(molecule("N", [(34, 0, 0)])))
    oxygen = next(a for a in peptide.GetAtoms() if a.GetSymbol() == "O" and any(b.GetBondTypeAsDouble() == 2 for b in a.GetBonds()))
    carbon = oxygen.GetNeighbors()[0]
    origin = np.array(peptide.GetConformer().GetAtomPosition(oxygen.GetIdx()))
    direction = np.array(peptide.GetConformer().GetAtomPosition(carbon.GetIdx())) - origin
    direction /= np.linalg.norm(direction)
    perpendicular = np.cross(direction, (0, 0, 1))
    perpendicular /= np.linalg.norm(perpendicular)
    axis = -.5 * direction + np.sqrt(.75) * perpendicular
    halogen = molecule("CCl", [tuple(origin+4.7*axis), tuple(origin+3*axis)])
    assert any(c["interaction_type"] == "XBDonor" and c["insertion_code"] == "A" for c in engine.analyze(halogen))
    # Reverse the C-Cl axis at fixed Cl position: geometric near miss.
    halogen.GetConformer().SetAtomPosition(0, tuple(origin+1.3*axis))
    assert not any(c["interaction_type"] == "XBDonor" for c in engine.analyze(halogen))


def test_postdock_optional_analysis_reuses_existing_exports(project, monkeypatch):
    from moldockpipe.pipeline import PipelineRunner
    repo, _ = project
    repo.update_settings({"postdock": {"mode": "sdf_only", "analyze_interactions": True}})
    calls = []
    runner = PipelineRunner(repo)
    monkeypatch.setattr(runner, "run_interactions", lambda **kw: calls.append(kw) or (2, 0))
    assert runner.run_postdock() == (0, 0)  # existing exports skipped
    assert calls == [{"profile_id": "r"}]
    repo.update_settings({"postdock": {"analyze_interactions": False}})
    runner.run_postdock()
    assert len(calls) == 1


def test_analysis_cli_never_imports_or_prepares(project, monkeypatch):
    from moldockpipe.__main__ import _run_project
    from moldockpipe.pipeline import PipelineRunner
    repo, _ = project
    def forbidden(*args, **kwargs):
        raise AssertionError("Upstream workflow invoked")
    monkeypatch.setattr(ProjectRepository, "import_project_input", forbidden)
    monkeypatch.setattr(PipelineRunner, "validate_docking_configuration", forbidden)
    monkeypatch.setattr(PipelineRunner, "run_interactions", lambda self, **kwargs: (2, 0))
    assert _run_project(repo.root, ["interactions"], profile_id="r") == 0


def test_calculation_failure_and_cancellation_mid_batch(project):
    repo, _ = project
    class PartialEngine(FakeEngine):
        def analyze(self, mol):
            if mol.GetConformer().GetAtomPosition(0).y:
                raise RuntimeError("candidate detector failure")
            return [CONTACT]
    assert InteractionAnalysisService(repo, engine_factory=PartialEngine).run(profile_id="r") == (1, 1)
    folder, rows = reports(repo)
    assert rows[1]["total_interaction_pairs"] == "NA"
    stop = [False]
    def progress(event):
        if event.event == "item_succeeded":
            stop[0] = True
    with pytest.raises(InterruptedError):
        service(repo, progress=progress, cancelled=lambda: stop[0]).run(profile_id="r", force=True)
    assert not list(folder.parent.glob(".pending-*"))
    with repo.connection() as conn:
        assert conn.execute("SELECT count(*) FROM interaction_runs WHERE status='cancelled'").fetchone()[0] == 1


def test_unreadable_file_does_not_abort_other_files(project, monkeypatch):
    import moldockpipe.services.interaction_analysis as module
    repo, path = project
    write_sdf(path.with_name("readable.sdf"), [molecule()])
    original = module.file_sha256
    def digest(p):
        if p.name == "all.sdf":
            raise PermissionError("fixture unreadable")
        return original(p)
    monkeypatch.setattr(module, "file_sha256", digest)
    assert service(repo).run(profile_id="r") == (1, 1)
    _, rows = reports(repo)
    assert "Unreadable SDF" in rows[0]["error_message"]


def test_missing_raw_poses_are_failed_records(project):
    repo, _ = project
    (repo.root / "artifacts/docking/raw.pdbqt").unlink()
    assert service(repo).run(profile_id="r") == (0, 2)
    _, rows = reports(repo)
    assert all("Source docking pose unavailable" in row["error_message"] for row in rows)


def test_corrupted_report_recomputed_and_failed_publication_preserves_latest(project, monkeypatch):
    repo, _ = project
    service(repo).run(profile_id="r")
    folder, _ = reports(repo)
    (folder / "analysis_manifest.json").write_text("corrupt")
    assert service(repo).run(profile_id="r") == (2, 0)
    current, _ = reports(repo)
    assert current != folder
    pointer = (folder.parent / "latest.json").read_bytes()
    original = Path.rename
    def rename(path, target):
        if path.name == "reports":
            raise OSError("publication failure")
        return original(path, target)
    monkeypatch.setattr(Path, "rename", rename)
    with pytest.raises(OSError, match="publication failure"):
        service(repo).run(profile_id="r", force=True)
    assert (folder.parent / "latest.json").read_bytes() == pointer
    assert not list(folder.parent.glob(".pending-*"))


def test_real_engine_full_job(project):
    from rdkit.Chem import AllChem
    repo, _ = project
    peptide = Chem.AddHs(Chem.MolFromSequence("A"), addResidueInfo=True)
    AllChem.EmbedMolecule(peptide, randomSeed=42)
    directory = repo.root / "inputs/receptors/r"
    protein = Chem.MolToPDBBlock(peptide).split("CONECT")[0].split("END")[0]
    for name in ("receptor_cleaned.pdb", "receptor_prepared.pdb", "receptor.pdbqt"):
        (directory / name).write_text(protein + "END\n")
    (directory / "source.pdb").write_text(protein + (directory / "reference.pdb").read_text())
    report = {key: file_sha256(directory / name) for key, name in (
        ("source_sha256", "source.pdb"), ("cleaned_receptor_sha256", "receptor_cleaned.pdb"),
        ("prepared_receptor_sha256", "receptor.pdbqt"))}
    (directory / "preparation_report.json").write_text(json.dumps(report))
    with repo.connection() as conn:
        conn.execute("UPDATE docking_runs SET receptor_hash=?", (report["prepared_receptor_sha256"],))
    assert InteractionAnalysisService(repo).run(profile_id="r") == (2, 0)
    _, rows = reports(repo)
    assert all(row["analysis_status"] == "success" for row in rows)
    assert int(rows[0]["total_interaction_pairs"]) > 0


def test_additive_migration_does_not_repeat_library_restoration(project, monkeypatch):
    repo, _ = project
    def unexpected(*args):
        raise AssertionError("Previous migration must not rerun")
    monkeypatch.setattr(repo, "_restore_archived_library", unexpected)
    repo.migrate()


def test_attach_crystal_reference_verifies_source_and_preserves_validation_metadata(project):
    from moldockpipe.services.interaction_analysis import attach_crystallographic_reference
    repo, _ = project
    before = repo.get_receptor_profiles()[0]["reference_ligand"]
    sdf = repo.root / before["sdf"]
    record = attach_crystallographic_reference(repo, "r", sdf)
    assert record["coordinate_validation"] == "source_coordinates_verified"
    assert repo.get_receptor_profiles()[0]["reference_ligand"] == before
    assert (repo.root / record["sdf"]).read_bytes() == sdf.read_bytes()
    assert service(repo).run(profile_id="r") == (2, 0)
    folder, _ = reports(repo)
    assert json.loads((folder / "analysis_manifest.json").read_text())["reference_origin"]["sha256"] == record["sha256"]


def test_imported_reference_without_validation_run_requires_declaration(project, tmp_path):
    from moldockpipe.services.interaction_analysis import attach_crystallographic_reference
    repo, _ = project
    external = tmp_path / "crystal.sdf"
    write_sdf(external, [molecule()])
    profiles = repo.get_receptor_profiles()
    profiles[0].pop("reference_ligand")
    repo.save_receptor_profiles(profiles)
    directory = repo.root / "inputs/receptors/r"
    for name in ("preparation_report.json", "source.pdb", "receptor_cleaned.pdb", "reference.pdb", "reference.sdf"):
        (directory / name).unlink()
    with pytest.raises(ValueError, match="Confirm"):
        attach_crystallographic_reference(repo, "r", external)
    record = attach_crystallographic_reference(repo, "r", external, confirm_same_frame=True)
    assert record["coordinate_validation"] == "user_confirmed_same_frame"
    external.unlink()  # the project owns an immutable copy
    assert service(repo).run(profile_id="r") == (2, 0)
    with repo.connection() as conn:
        assert conn.execute("SELECT count(*) FROM redocking_runs").fetchone()[0] == 0
        assert conn.execute("SELECT count(*) FROM docking_runs WHERE is_current=1").fetchone()[0] == 1
    (repo.root / record["sdf"]).write_text("changed")
    with pytest.raises(ValueError, match="has changed"):
        service(repo).run(profile_id="r")


def test_manual_reference_rejects_wrong_frame_invalid_and_multiple_records(project, tmp_path):
    from moldockpipe.services.interaction_analysis import attach_crystallographic_reference
    repo, _ = project
    before = repo.get_receptor_profiles()
    path = tmp_path / "reference.sdf"
    write_sdf(path, [molecule(points=[(100, 0, 0), (101.4, 0, 0)])])
    with pytest.raises(ValueError, match="source coordinates"):
        attach_crystallographic_reference(repo, "r", path, confirm_same_frame=True)
    write_sdf(path, [molecule(), molecule()])
    with pytest.raises(ValueError, match="exactly one"):
        attach_crystallographic_reference(repo, "r", path)
    path.write_text("invalid\n$$$$\n")
    with pytest.raises(ValueError, match="Invalid"):
        attach_crystallographic_reference(repo, "r", path)
    assert repo.get_receptor_profiles() == before


def test_manual_reference_verifies_receptor_frame_too(project):
    from moldockpipe.services.interaction_analysis import attach_crystallographic_reference
    repo, _ = project
    directory = repo.root / "inputs/receptors/r"
    for name in ("receptor.pdbqt", "receptor_prepared.pdb"):
        path = directory / name
        path.write_text(path.read_text().replace("10.000", "20.000"))
    report_path = directory / "preparation_report.json"
    report = json.loads(report_path.read_text())
    report["prepared_receptor_sha256"] = file_sha256(directory / "receptor.pdbqt")
    report_path.write_text(json.dumps(report))
    with pytest.raises(ValueError, match="Prepared receptor coordinates differ"):
        attach_crystallographic_reference(repo, "r", directory / "reference.sdf", confirm_same_frame=True)
