from pathlib import Path
import hashlib

import pytest

from moldockpipe.project import ProjectRepository, utc_now
from moldockpipe.services.structure_viewer import StructureViewerService, extract_pose


ATOM = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00"
PDB = ATOM.ljust(76) + " C\nEND\n"
PDBQT = ATOM.ljust(70) + "-0.123 C\n"
SDF = "pose\n  test\n\n  1  0  0  0  0  0            999 V2000\n    1.0000    2.0000    3.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n"


@pytest.fixture
def viewer_project(tmp_path, monkeypatch):
    repo = ProjectRepository.create(tmp_path / "project")
    csv = tmp_path / "ligands.csv"
    csv.write_text("id,smiles\nligand,CC\n", encoding="utf-8")
    repo.import_ligands_csv(csv)
    folder = repo.root / "inputs/receptors/target"
    folder.mkdir()
    receptor = folder / "receptor.pdbqt"
    receptor.write_text(PDBQT, encoding="utf-8")
    (folder / "receptor_prepared.pdb").write_text(PDB, encoding="utf-8")
    repo.save_receptor_profiles([{"id": "target", "name": "Target", "enabled": True,
                                  "receptor": receptor.relative_to(repo.root).as_posix()}])
    raw = repo.root / "artifacts/docking/raw.pdbqt"
    raw.write_text("REMARK SMILES CC\nMODEL 3\nREMARK SMILES IDX 1 1\n" + PDBQT +
                   "ENDMDL\nMODEL 8\nREMARK SMILES IDX 1 1\n" + PDBQT.replace("1.000", "8.000") + "ENDMDL\n", encoding="utf-8")
    artifact = repo.add_artifact(raw, "vina_output_pdbqt", "vina")
    with repo.connection() as conn:
        conn.execute("""INSERT INTO molecular_states
            (state_id,parent_id,state_smiles,state_isomeric_smiles,formal_charge,state_structure_hash,generation_fingerprint,status,created_at)
            VALUES ('state','ligand','CC','CC',0,'structure','generation','prepared',?)""", (utc_now(),))
        conn.execute("""INSERT INTO docking_runs
            (run_id,state_id,receptor_hash,settings_fingerprint,status,raw_output_artifact_id,
             is_current,receptor_profile_id,receptor_profile_name,command_json)
            VALUES ('run','state',?,'settings','completed',?,1,'target','Target','{}')""",
                     (hashlib.sha256(receptor.read_bytes()).hexdigest(), artifact))
        for mode, score in ((3, -5), (8, -9)):
            conn.execute("INSERT INTO docking_poses(pose_id,run_id,mode_index,affinity) VALUES (?,'run',?,?)",
                         (f"pose-{mode}", mode, score))
    service = StructureViewerService(repo)
    calls = []

    def convert(source, output, log, **kwargs):
        calls.append((source, source.read_text(encoding="utf-8"), output.parent))
        output.write_text(SDF, encoding="utf-8")

    monkeypatch.setattr(service.postdock, "export_sdf", convert)
    return repo, service, calls, receptor, raw


def test_extracts_actual_mode_with_global_and_local_meeko_remarks():
    raw = "REMARK SMILES CC\nMODEL 8\nREMARK SMILES IDX 1 1\n" + PDBQT + "ENDMDL\nMODEL 3\n" + PDBQT.replace("1.000", "3.000") + "ENDMDL\n"
    selected = extract_pose(raw, 8)
    assert selected.startswith("REMARK SMILES CC\nMODEL 8\nREMARK SMILES IDX 1 1\n")
    assert PDBQT in selected
    assert "3.000   2.000" not in selected
    with pytest.raises(ValueError, match="mode 1"):
        extract_pose(raw, 1)


def test_rejects_ambiguous_models_and_handles_single_pose():
    with pytest.raises(ValueError, match="one valid"):
        extract_pose(f"MODEL 1\n{PDBQT}ENDMDL\nMODEL 1\n{PDBQT}ENDMDL\n", 1)
    assert extract_pose(PDBQT, 1) == PDBQT
    with pytest.raises(ValueError):
        extract_pose(PDBQT, 2)


def test_load_uses_exact_pose_preserves_coordinates_caches_and_cleans(viewer_project):
    repo, service, calls, receptor, raw = viewer_project
    result = service.load("run", "pose-8")
    assert result.mode_index == 8 and result.affinity == -9
    assert result.run_id == "run" and result.state_id == "state"
    assert result.receptor_pdb == PDB
    assert "8.000   2.000   3.000" in calls[0][1]
    assert "1.000   2.000   3.000" not in calls[0][1]
    assert not calls[0][0].exists() and not calls[0][2].exists()
    assert service.load("run", "pose-8") == result
    assert len(calls) == 1
    assert [p["mode_index"] for p in service.poses("run")] == [3, 8]
    service.load("run", "pose-3")
    assert len(calls) == 2


@pytest.mark.parametrize("which", ["receptor", "raw", "prepared"])
def test_rejects_changed_or_corrupt_inputs_even_when_cached(viewer_project, which):
    repo, service, calls, receptor, raw = viewer_project
    service.load("run", "pose-8")
    path = {"receptor": receptor, "raw": raw, "prepared": receptor.parent / "receptor_prepared.pdb"}[which]
    path.write_text(path.read_text().replace("1.000", "9.000"), encoding="utf-8")
    with pytest.raises(ValueError, match="changed|checksum|does not match"):
        service.load("run", "pose-8")
    assert len(calls) == 1


def test_imported_receptor_falls_back_to_pdbqt_never_cleaned(viewer_project):
    repo, service, calls, receptor, raw = viewer_project
    (receptor.parent / "receptor_prepared.pdb").unlink()
    (receptor.parent / "receptor_cleaned.pdb").write_text(PDB.replace("1.000", "9.000"))
    result = service.load("run", "pose-3")
    assert "1.000   2.000   3.000" in result.receptor_pdb
    assert "9.000" not in result.receptor_pdb


def test_load_finds_exact_archived_artifact_and_verifies_checksum(viewer_project):
    repo, service, calls, receptor, raw = viewer_project
    relative = raw.relative_to(repo.root)
    archived = repo.root / "Run3" / relative
    archived.parent.mkdir(parents=True)
    raw.rename(archived)
    wrong = repo.root / "Run1" / relative
    wrong.parent.mkdir(parents=True)
    wrong.write_text("different docking run", encoding="utf-8")
    assert service.load("run", "pose-8").mode_index == 8
    assert "8.000   2.000   3.000" in calls[0][1]
    assert not raw.exists()  # Viewing never moves files or rewrites the record.
    archived.write_text("corrupt archive", encoding="utf-8")
    with pytest.raises(ValueError, match="checksums"):
        service.load("run", "pose-8")


def test_corrupt_recorded_artifact_cannot_use_good_archive(viewer_project):
    repo, service, calls, receptor, raw = viewer_project
    archived = repo.root / "Run3" / raw.relative_to(repo.root)
    archived.parent.mkdir(parents=True)
    archived.write_bytes(raw.read_bytes())
    raw.write_text("corrupt current file", encoding="utf-8")
    with pytest.raises(ValueError, match="checksum"):
        service.load("run", "pose-8")
    assert not calls


def test_archive_lookup_requires_complete_relative_path(viewer_project):
    repo, service, calls, receptor, raw = viewer_project
    wrong = repo.root / "Run3" / "another-run" / raw.name
    wrong.parent.mkdir(parents=True)
    raw.rename(wrong)
    with pytest.raises(FileNotFoundError, match="Restore the original artifact"):
        service.load("run", "pose-8")
    assert not calls


@pytest.mark.parametrize("relative", ["../outside.pdbqt", "artifacts/../../outside.pdbqt"])
def test_artifact_lookup_rejects_path_traversal(viewer_project, relative):
    repo, service, calls, receptor, raw = viewer_project
    with pytest.raises(ValueError, match="outside this project"):
        service._raw_artifact(relative, hashlib.sha256(raw.read_bytes()).hexdigest())


def test_retained_metal_is_preserved_without_scoring_sites(viewer_project):
    repo, service, calls, receptor, raw = viewer_project
    metal = "HETATM    2 ZN    ZN A   2       4.000   5.000   6.000  1.00 20.00"
    receptor.write_text(PDBQT + metal.ljust(76) + "Zn\n")
    (receptor.parent / "receptor_prepared.pdb").write_text(PDB + metal.ljust(76) + "Zn\n")
    with repo.connection() as conn:
        conn.execute("UPDATE docking_runs SET receptor_hash=?", (hashlib.sha256(receptor.read_bytes()).hexdigest(),))
    assert "Zn" in service.load("run", "pose-8").receptor_pdb


def test_conversion_failure_cleans_temporary_files(viewer_project, monkeypatch):
    repo, service, calls, receptor, raw = viewer_project
    directories = []

    def fail(source, output, log):
        directories.append(source.parent)
        raise RuntimeError("conversion failed")

    monkeypatch.setattr(service.postdock, "export_sdf", fail)
    with pytest.raises(RuntimeError, match="conversion failed"):
        service.load("run", "pose-8")
    assert not directories[0].exists()
    assert not service._cache


def test_interaction_analysis_is_cached_and_runs_with_original_coordinates(viewer_project, monkeypatch):
    from moldockpipe.services import interactions
    repo, service, calls, receptor, raw = viewer_project
    analyses = []

    def analyze(pdb, sdf, *, cancelled):
        analyses.append((pdb, sdf))
        return [], "Geometric predictions"

    monkeypatch.setattr(interactions, "analyze_interactions", analyze)
    result = service.load("run", "pose-8")
    assert result.interaction_notes == "Geometric predictions"
    assert analyses == [(PDB, SDF)]
    service.load("run", "pose-8")
    assert len(analyses) == 1
    service.load("run", "pose-3")
    assert len(analyses) == 2


def test_interaction_failure_keeps_structure_viewable(viewer_project, monkeypatch):
    from moldockpipe.services import interactions
    repo, service, calls, receptor, raw = viewer_project

    def fail(*args, **kwargs):
        raise ValueError("chemistry unavailable")

    monkeypatch.setattr(interactions, "analyze_interactions", fail)
    result = service.load("run", "pose-8")
    assert result.receptor_pdb == PDB and result.ligand_sdf == SDF
    assert result.interactions == []
    assert "Interaction analysis unavailable: chemistry unavailable" == result.interaction_notes


def test_missing_files_and_cancelled_load(viewer_project):
    repo, service, calls, receptor, raw = viewer_project
    with pytest.raises(InterruptedError):
        service.load("run", "pose-8", cancelled=lambda: True)
    raw.unlink()
    with pytest.raises(FileNotFoundError):
        service.load("run", "pose-8")
    with pytest.raises(ValueError):
        service.load("run", "other-pose")


def test_large_receptor_is_not_embedded_in_html(viewer_project):
    repo, service, calls, receptor, raw = viewer_project
    # A >2MB receptor is handed across the bridge, never encoded into setHtml().
    count = 30000
    receptor.write_text(PDBQT * count)
    (receptor.parent / "receptor_prepared.pdb").write_text(PDB * count)
    with repo.connection() as conn:
        conn.execute("UPDATE docking_runs SET receptor_hash=?", (hashlib.sha256(receptor.read_bytes()).hexdigest(),))
    assert len(service.load("run", "pose-8").receptor_pdb) > 2 * 1024 * 1024


def test_real_meeko_conversion_keeps_selected_heavy_atom_coordinates(tmp_path):
    meeko = pytest.importorskip("meeko")
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from moldockpipe.services.postdock import PostDockService, _pdbqt_element
    molecule = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    assert AllChem.EmbedMolecule(molecule, randomSeed=42) == 0
    setup = meeko.MoleculePreparation().prepare(molecule)[0]
    pdbqt, success, error = meeko.PDBQTWriterLegacy.write_string(setup)
    assert success, error
    raw = f"MODEL 2\n{pdbqt}ENDMDL\nMODEL 8\n{pdbqt}ENDMDL\n"
    source = tmp_path / "pose.pdbqt"
    source.write_text(extract_pose(raw, 8), encoding="utf-8")
    output = tmp_path / "pose.sdf"
    PostDockService().export_sdf(source, output, tmp_path / "conversion.log", cancelled=lambda: False)
    molecules = [mol for mol in Chem.SDMolSupplier(str(output), removeHs=False) if mol is not None]
    assert len(molecules) == 1
    actual = sorted(tuple(molecules[0].GetConformer().GetAtomPosition(atom.GetIdx()))
                    for atom in molecules[0].GetAtoms() if atom.GetAtomicNum() > 1)
    expected = sorted(tuple(float(line[start:start + 8]) for start in (30, 38, 46))
                      for line in pdbqt.splitlines()
                      if line.startswith(("ATOM  ", "HETATM")) and _pdbqt_element(line) != "H")
    assert len(actual) == len(expected)
    for displayed, docked in zip(actual, expected):
        assert displayed == pytest.approx(docked, abs=0.0001)


@pytest.mark.parametrize("failure", ["cancel", "timeout"])
def test_viewer_conversion_terminates_child_and_cleans_up(tmp_path, monkeypatch, failure):
    import moldockpipe.services.postdock as module
    source = tmp_path / "pose.pdbqt"
    source.write_text(PDBQT)
    processes = []

    class Process:
        def __init__(self, command, *, cwd, **kwargs):
            self.cwd, self.returncode, self.killed = cwd, None, False
            processes.append(self)

        def poll(self):
            return self.returncode

        def kill(self):
            self.killed = True
            self.returncode = -1

        def communicate(self, **kwargs):
            return "", ""

    monkeypatch.setattr(module.subprocess, "Popen", Process)
    if failure == "timeout":
        clock = iter((0, 121))
        monkeypatch.setattr(module.time, "monotonic", lambda: next(clock))
    with pytest.raises(InterruptedError if failure == "cancel" else TimeoutError):
        module.PostDockService().export_sdf(source, tmp_path / "out.sdf", tmp_path / "out.log",
                                           cancelled=lambda: failure == "cancel")
    assert processes[0].killed and not processes[0].cwd.exists()
