from pathlib import Path

import pytest

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


def test_ligand_sync_is_idempotent_and_incremental(tmp_path: Path) -> None:
    repo = ProjectRepository.create(tmp_path / "example")
    source = tmp_path / "ligands.csv"
    source.write_text("id,smiles,notes,params_json\na,CCO,first,{}\nb,CCN,second,{}\n", encoding="utf-8")

    first = repo.sync_ligands_csv(source)
    second = repo.sync_ligands_csv(source)
    assert (first.added, first.unchanged, first.changed, first.archived) == (2, 0, 0, 0)
    assert (second.added, second.unchanged, second.changed, second.archived) == (0, 2, 0, 0)

    with repo.connection() as conn:
        conn.execute("""INSERT INTO screening_results
            (parent_id, fingerprint, status, descriptors_json, rules_json, decision, reason, created_at, active)
            VALUES ('a', 'old', 'completed', '{}', '{}', 'warning', NULL, 'now', 1)""")

    source.write_text("id,smiles,notes,params_json\na,CCC,changed,{}\nc,COC,new,{}\n", encoding="utf-8")
    result = repo.sync_ligands_csv(source)
    assert (result.added, result.unchanged, result.changed, result.archived) == (1, 0, 1, 1)
    with repo.connection() as conn:
        rows = conn.execute("SELECT parent_id, active, source_smiles FROM parent_ligands ORDER BY parent_id").fetchall()
        screening = conn.execute("SELECT active FROM screening_results WHERE parent_id='a'").fetchone()[0]
    assert [tuple(row) for row in rows] == [("a", 1, "CCC"), ("b", 0, "CCN"), ("c", 1, "COC")]
    assert screening == 0


def test_ligand_sync_rejects_missing_or_duplicate_ids_atomically(tmp_path: Path) -> None:
    repo = ProjectRepository.create(tmp_path / "example")
    valid = tmp_path / "valid.csv"
    valid.write_text("id,smiles\na,CCO\n", encoding="utf-8")
    repo.sync_ligands_csv(valid)

    missing_id = tmp_path / "missing.csv"
    missing_id.write_text("id,smiles\n,CCN\n", encoding="utf-8")
    with pytest.raises(ValueError, match="stable ligand id"):
        repo.sync_ligands_csv(missing_id)

    duplicate_id = tmp_path / "duplicate.csv"
    duplicate_id.write_text("id,smiles\na,CCN\na,CCC\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Duplicate ligand id"):
        repo.sync_ligands_csv(duplicate_id)

    with repo.connection() as conn:
        assert conn.execute("SELECT source_smiles FROM parent_ligands WHERE parent_id='a'").fetchone()[0] == "CCO"
