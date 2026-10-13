import sqlite3

import pytest

from moldockpipe.pipeline import PipelineRunner
from moldockpipe.project import ProjectRepository


def test_new_project_includes_empty_input_csv(tmp_path) -> None:
    root = tmp_path / "project"

    ProjectRepository.create(root)

    input_csv = root / "inputs" / "input.csv"
    assert input_csv.is_file()
    assert input_csv.read_text(encoding="utf-8") == "id,smiles,notes,params_json\n"


def test_project_creation_preserves_existing_input_csv(tmp_path) -> None:
    root = tmp_path / "project"
    input_csv = root / "inputs" / "input.csv"
    input_csv.parent.mkdir(parents=True)
    input_csv.write_text("id,smiles,notes,params_json\nlig1,CCO,,\n", encoding="utf-8")

    ProjectRepository.create(root)

    assert input_csv.read_text(encoding="utf-8") == "id,smiles,notes,params_json\nlig1,CCO,,\n"


def test_repeated_chemistry_stages_reuse_artifacts(tmp_path, monkeypatch) -> None:
    pytest.importorskip("rdkit")
    pytest.importorskip("molscrub")
    pytest.importorskip("meeko")
    repo = ProjectRepository.create(tmp_path / "project")
    csv_path = tmp_path / "input.csv"
    csv_path.write_text("id,smiles\nlig1,CCO\n", encoding="utf-8")
    repo.import_ligands_csv(csv_path)
    runner = PipelineRunner(repo)
    runner.run_screening()
    runner.run_molscrub()
    runner.run_meeko()
    # Reimporting and editing annotations must not call chemistry tools again.
    csv_path.write_text("id,smiles,notes\nlig1,CCO,updated note\n", encoding="utf-8")
    assert repo.sync_ligands_csv(csv_path).changed == 1
    def unexpected_work(*args, **kwargs):
        pytest.fail("Completed chemistry was processed again")
    with monkeypatch.context() as cached:
        cached.setattr("moldockpipe.pipeline.screen_smiles", unexpected_work)
        cached.setattr("moldockpipe.pipeline.MolScrubService.generate", unexpected_work)
        cached.setattr("moldockpipe.pipeline.MeekoService.prepare_ligand", unexpected_work)
        runner.run_screening()
        runner.run_molscrub()
        runner.run_meeko()
    with sqlite3.connect(repo.database_path) as conn:
        assert conn.execute("SELECT COUNT(*) FROM parent_ligands").fetchone()[0] == 1
        assert conn.execute("SELECT COUNT(*) FROM molecular_states").fetchone()[0] == 1
        assert conn.execute("SELECT COUNT(*) FROM artifacts").fetchone()[0] == 3


@pytest.mark.parametrize("damage", ["missing", "corrupt", "failed"])
def test_prepared_file_is_repaired_without_regenerating_states(tmp_path, monkeypatch, damage) -> None:
    pytest.importorskip("rdkit")
    pytest.importorskip("molscrub")
    pytest.importorskip("meeko")
    repo = ProjectRepository.create(tmp_path / "project")
    source = repo.root / "inputs" / "input.csv"
    source.write_text("id,smiles\nlig1,CCO\n", encoding="utf-8")
    repo.import_project_input()
    runner = PipelineRunner(repo)
    runner.run_screening()
    runner.run_molscrub()
    assert runner.run_meeko() == (1, 0)
    with repo.connection() as conn:
        path = conn.execute("SELECT relative_path FROM artifacts WHERE artifact_type='ligand_pdbqt'").fetchone()[0]
    if damage == "missing":
        (repo.root / path).unlink()
    elif damage == "corrupt":
        (repo.root / path).write_text("damaged", encoding="utf-8")
    else:
        with repo.connection() as conn:
            conn.execute("UPDATE conformers SET status='pdbqt_failed'")
    def unexpected_generation(*args, **kwargs):
        pytest.fail("Repairing a PDBQT should not regenerate valid SDF states")
    monkeypatch.setattr("moldockpipe.pipeline.MolScrubService.generate", unexpected_generation)
    assert runner.run_molscrub() == (1, 0)
    assert runner.run_meeko() == (1, 0)
    assert (repo.root / path).is_file()


def test_state_regeneration_preserves_old_docking_history(tmp_path) -> None:
    pytest.importorskip("rdkit")
    pytest.importorskip("molscrub")
    repo = ProjectRepository.create(tmp_path / "project")
    source = repo.root / "inputs" / "input.csv"
    source.write_text("id,smiles\nlig1,CCO\n", encoding="utf-8")
    repo.import_project_input()
    runner = PipelineRunner(repo)
    runner.run_screening()
    runner.run_molscrub()
    with repo.connection() as conn:
        state_id = conn.execute("SELECT state_id FROM molecular_states").fetchone()[0]
        conn.execute("""INSERT INTO docking_runs (run_id,state_id,receptor_hash,settings_fingerprint,status,is_current)
            VALUES ('old-run',?,'receptor','settings','completed',1)""", (state_id,))
        conn.execute("INSERT INTO docking_poses (pose_id,run_id,mode_index,affinity) VALUES ('pose','old-run',1,-7)")
        path = conn.execute("SELECT relative_path FROM artifacts WHERE artifact_type='state_sdf'").fetchone()[0]
    (repo.root / path).unlink()
    assert runner.run_molscrub() == (1, 0)
    with repo.connection() as conn:
        assert conn.execute("SELECT COUNT(*) FROM docking_poses").fetchone()[0] == 1
        assert conn.execute("SELECT is_current FROM docking_runs WHERE run_id='old-run'").fetchone()[0] == 0
        assert conn.execute("PRAGMA foreign_key_check").fetchall() == []
