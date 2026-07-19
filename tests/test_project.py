from pathlib import Path

import pytest

from moldockpipe.project import ProjectRepository


def test_project_creation_and_csv_import(tmp_path: Path) -> None:
    repo = ProjectRepository.create(tmp_path / "example")
    assert repo.get_receptor_profiles() == []
    assert not (repo.root / "inputs" / "receptors" / "default").exists()
    source = tmp_path / "ligands.csv"
    source.write_text("id,smiles,notes,params_json\nlig1,CCO,ethanol,{}\n", encoding="utf-8")
    assert repo.import_ligands_csv(source) == 1
    with repo.connection() as conn:
        row = conn.execute("SELECT source_smiles, parse_status FROM parent_ligands WHERE parent_id='lig1'").fetchone()
    assert tuple(row) == ("CCO", "pending")


def test_remove_receptor_profile_deletes_only_its_receptor_folder(tmp_path: Path) -> None:
    repo = ProjectRepository.create(tmp_path / "example")
    receptor_folder = repo.root / "inputs" / "receptors" / "to-delete"
    receptor_folder.mkdir(parents=True)
    (receptor_folder / "receptor.pdbqt").write_text("RECEPTOR", encoding="utf-8")
    repo.save_receptor_profiles([{"id": "to-delete", "name": "Delete me", "enabled": True, "archived": False}])
    assert repo.remove_receptor_profile("to-delete") is True
    assert repo.get_receptor_profiles() == []
    assert not receptor_folder.exists()
    assert (repo.root / "inputs" / "receptors").exists()


def test_recover_interrupted_stage_runs(tmp_path: Path) -> None:
    repo = ProjectRepository.create(tmp_path / "example")
    with repo.connection() as conn:
        conn.execute("INSERT INTO stage_runs VALUES ('run', 'screening', 'parent', 'p1', 'fp', 'running', 'start', NULL, NULL)")
    assert repo.recover_interrupted_runs() == 1
    with repo.connection() as conn:
        assert conn.execute("SELECT status FROM stage_runs WHERE stage_run_id='run'").fetchone()[0] == "interrupted"


def test_redocking_validation_queue_is_persistent_fifo_and_resumable(tmp_path: Path) -> None:
    repo = ProjectRepository.create(tmp_path / "example")
    first = {"id": "rec-a", "name": "Receptor A"}
    second = {"id": "rec-b", "name": "Receptor B"}
    settings = {"exhaustiveness": 32, "num_modes": 20, "energy_range": 5.0, "seed": 123, "cpu_count": 2}
    first_id = repo.enqueue_redocking(first, settings)
    second_id = repo.enqueue_redocking(second, settings)
    with pytest.raises(ValueError, match="already has"):
        repo.enqueue_redocking(first, settings)

    claimed = repo.claim_next_redocking()
    assert claimed and claimed["queue_id"] == first_id and claimed["status"] == "running"
    assert claimed["redocking_run_id"]
    repo.update_redocking_queue_stage(first_id, "vina_redocking")

    # Opening the project after an interruption places the active item back at
    # the front with the same run ID, so completed redocking stages can be reused.
    assert repo.recover_interrupted_runs() >= 1
    resumed = repo.claim_next_redocking()
    assert resumed and resumed["queue_id"] == first_id
    assert resumed["redocking_run_id"] == claimed["redocking_run_id"]
    repo.finish_redocking_queue_item(first_id, "completed")
    next_item = repo.claim_next_redocking()
    assert next_item and next_item["queue_id"] == second_id


def test_redocking_queue_cancel_retry_and_clear(tmp_path: Path) -> None:
    repo = ProjectRepository.create(tmp_path / "example")
    profile = {"id": "rec-a", "name": "Receptor A"}
    settings = {"exhaustiveness": 32, "num_modes": 20, "energy_range": 5.0, "seed": 123, "cpu_count": 2}
    queue_id = repo.enqueue_redocking(profile, settings)
    assert repo.cancel_redocking_queue_item(queue_id)
    assert repo.list_redocking_queue()[0]["status"] == "cancelled"
    assert repo.retry_redocking_queue_item(queue_id)
    assert repo.claim_next_redocking()["queue_id"] == queue_id
    repo.finish_redocking_queue_item(queue_id, "completed")
    assert repo.clear_finished_redocking_queue() == 1
    assert repo.list_redocking_queue() == []


def test_structured_workflow_and_receptor_provenance(tmp_path: Path) -> None:
    repo = ProjectRepository.create(tmp_path / "example")
    workflow_id = repo.start_workflow_run("docking", {"receptors": ["rec-a"], "requested": 10})
    repo.update_workflow_stage(workflow_id, "vina", "running", summary={"total": 10})
    repo.record_provenance_event(workflow_run_id=workflow_id, event_type="campaign_started",
                                 stage_name="vina", data={"total": 10})
    repo.update_workflow_stage(workflow_id, "vina", "completed", summary={"succeeded": 9, "failed": 1})
    repo.finish_workflow_run(workflow_id, "completed")

    source = tmp_path / "source.pdb"; source.write_text("ATOM", encoding="ascii")
    prep_id = repo.start_receptor_preparation(profile_id="rec-a", receptor_name="Receptor A",
        source_path=source, source_sha256="source-hash", inventory={"protein_residues": 100},
        decisions={"center": {"method": "reference_ligand_centroid"}})
    repo.finish_receptor_preparation(prep_id, "completed", command={"argv": ["meeko"]},
                                     artifacts={"receptor.pdbqt": "hash"}, warnings=["review altloc"])
    with repo.connection() as conn:
        workflow = conn.execute("SELECT status,settings_json FROM workflow_runs WHERE workflow_run_id=?", (workflow_id,)).fetchone()
        stage = conn.execute("SELECT status,summary_json FROM workflow_run_stages WHERE workflow_run_id=?", (workflow_id,)).fetchone()
        prep = conn.execute("SELECT status,inventory_json,decisions_json FROM receptor_preparation_runs WHERE preparation_run_id=?", (prep_id,)).fetchone()
        events = conn.execute("SELECT COUNT(*) FROM provenance_events WHERE workflow_run_id=?", (workflow_id,)).fetchone()[0]
    assert workflow["status"] == "completed" and '"requested": 10' in workflow["settings_json"]
    assert stage["status"] == "completed" and '"succeeded": 9' in stage["summary_json"]
    assert prep["status"] == "completed" and '"protein_residues": 100' in prep["inventory_json"]
    assert "reference_ligand_centroid" in prep["decisions_json"] and events == 1


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
