from pathlib import Path

from moldockpipe.project import ProjectRepository


def test_project_creation_and_csv_import(tmp_path: Path) -> None:
    repo = ProjectRepository.create(tmp_path / "example")
    source = tmp_path / "ligands.csv"
    source.write_text("id,smiles,notes,params_json\nlig1,CCO,ethanol,{}\n", encoding="utf-8")
    assert repo.import_ligands_csv(source) == 1
    with repo.connection() as conn:
        row = conn.execute("SELECT source_smiles, parse_status FROM parent_ligands WHERE parent_id='lig1'").fetchone()
    assert tuple(row) == ("CCO", "pending")


def test_recover_interrupted_stage_runs(tmp_path: Path) -> None:
    repo = ProjectRepository.create(tmp_path / "example")
    with repo.connection() as conn:
        conn.execute("INSERT INTO stage_runs VALUES ('run', 'screening', 'parent', 'p1', 'fp', 'running', 'start', NULL, NULL)")
    assert repo.recover_interrupted_runs() == 1
    with repo.connection() as conn:
        assert conn.execute("SELECT status FROM stage_runs WHERE stage_run_id='run'").fetchone()[0] == "interrupted"


def test_clear_workflow_removes_old_ligand_ids(tmp_path: Path) -> None:
    repo = ProjectRepository.create(tmp_path / "example")
    source = tmp_path / "ligands.csv"
    source.write_text("id,smiles\nold,CCO\n", encoding="utf-8")
    repo.import_ligands_csv(source)
    repo.clear_workflow_data()
    with repo.connection() as conn:
        assert conn.execute("SELECT COUNT(*) FROM parent_ligands").fetchone()[0] == 0
