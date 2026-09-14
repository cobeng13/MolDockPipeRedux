from __future__ import annotations

import csv
from pathlib import Path

import yaml

from moldockpipe.pipeline import PipelineRunner
from moldockpipe.project import ProjectRepository, utc_now
from moldockpipe.csv_exports import export_leaderboard_csv, export_manifest_csv


def test_legacy_vina_settings_migrate_to_default_profile(tmp_path: Path) -> None:
    root = tmp_path / "legacy"
    (root / "inputs").mkdir(parents=True)
    (root / "inputs" / "receptor_prepared.pdbqt").write_text("RECEPTOR", encoding="utf-8")
    (root / "project.yml").write_text(yaml.safe_dump({
        "schema_version": 2,
        "vina": {"receptor": "inputs/receptor_prepared.pdbqt", "center_x": 7.5, "size_x": 22, "cpu_count": 3},
    }), encoding="utf-8")

    repo = ProjectRepository(root)
    repo.migrate()

    profiles = repo.get_receptor_profiles()
    assert len(profiles) == 1
    assert profiles[0]["id"] == "default"
    assert profiles[0]["center_x"] == 7.5
    assert profiles[0]["cpu_count"] == 3
    assert (root / profiles[0]["receptor"]).read_text(encoding="utf-8") == "RECEPTOR"
    assert repo.get_settings()["schema_version"] == 6


def test_multidock_runs_and_reuses_each_receptor_independently(tmp_path: Path, monkeypatch) -> None:
    repo = ProjectRepository.create(tmp_path / "project")
    csv_path = tmp_path / "ligands.csv"
    csv_path.write_text("id,smiles\nlig1,CCO\n", encoding="utf-8")
    repo.import_ligands_csv(csv_path)

    state_id = "state-1"
    sdf = repo.root / "artifacts" / "sdf" / "lig1" / state_id / "state.sdf"
    sdf.parent.mkdir(parents=True, exist_ok=True); sdf.write_text("SDF", encoding="utf-8")
    sdf_artifact = repo.add_artifact(sdf, "state_sdf", "molscrub")
    ligand = repo.root / "artifacts" / "pdbqt" / "lig1" / state_id / "ligand.pdbqt"
    ligand.parent.mkdir(parents=True, exist_ok=True); ligand.write_text("LIGAND", encoding="utf-8")
    with repo.connection() as conn:
        conn.execute("""INSERT INTO molecular_states
            (state_id,parent_id,state_smiles,state_isomeric_smiles,state_inchikey,formal_charge,tautomer_index,protomer_index,
             state_structure_hash,generation_fingerprint,status,reason,enumeration_truncated,created_at,active)
            VALUES (?, 'lig1', 'CCO', 'CCO', NULL, 0, 0, 0, 'structure', 'generation', 'prepared', NULL, 0, ?, 1)""", (state_id, utc_now()))
        conn.execute("INSERT INTO conformers(conformer_id,state_id,conformer_index,sdf_artifact_id,status,created_at) VALUES ('conf-1', ?, 0, ?, 'pdbqt_ready', ?)", (state_id, sdf_artifact, utc_now()))

    profiles = []
    for profile_id, score in (("rec-a", -7.1), ("rec-b", -8.2)):
        receptor = repo.root / "inputs" / "receptors" / profile_id / "receptor.pdbqt"
        receptor.parent.mkdir(parents=True, exist_ok=True); receptor.write_text(profile_id, encoding="utf-8")
        (receptor.parent / "receptor_prepared.pdb").write_text(
            "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nEND\n",
            encoding="ascii",
        )
        profiles.append({
            "id": profile_id, "name": profile_id.upper(), "enabled": True, "archived": False,
            "receptor": receptor.relative_to(repo.root).as_posix(),
            "center_x": 0, "center_y": 0, "center_z": 0, "size_x": 20, "size_y": 20, "size_z": 20,
            "exhaustiveness": 8, "num_modes": 9, "energy_range": 3, "seed": 42, "cpu_count": 1,
            "test_score": score,
        })
    repo.save_receptor_profiles(profiles)
    standard_mol2 = """@<TRIPOS>MOLECULE
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
    with repo.connection() as conn:
        for profile in profiles:
            validation_run = f"validation-{profile['id']}"
            standard = repo.root / "inputs" / "receptors" / profile["id"] / "redocking" / validation_run / "dockrmsd" / "reference_ligand_heavy.mol2"
            standard.parent.mkdir(parents=True); standard.write_text(standard_mol2, encoding="utf-8")
            conn.execute("""INSERT INTO redocking_runs
                (run_id,receptor_profile_id,status,reference_ligand_id,receptor_path,receptor_sha256,
                 reference_sdf_path,reference_mol2_path,reference_sha256,settings_json,finished_at)
                VALUES (?, ?, 'ARTIFACTS_READY', 'LIG', 'receptor.pdbqt', 'receptor-hash',
                        'reference.sdf', 'reference.mol2', 'reference-hash', '{}', ?)""",
                (validation_run, profile["id"], utc_now()))

    calls: list[str] = []

    class FakeResult:
        def __init__(self, score: float) -> None:
            self.return_code = 0
            self.poses = [(1, score, 0.0, 0.0)]

    class FakeBackend:
        def run(self, *, receptor, output, log, settings, **kwargs):
            profile_id = receptor.parent.name
            calls.append(profile_id)
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text("POSE", encoding="utf-8"); log.write_text("LOG", encoding="utf-8")
            return FakeResult(float(settings["test_score"]))

    monkeypatch.setattr("moldockpipe.pipeline.VinaDockingBackend", FakeBackend)
    events = []
    runner = PipelineRunner(repo, progress=events.append)
    assert runner.run_vina() == (2, 0)
    assert calls == ["rec-a", "rec-b"]
    assert {event.receptor_profile_id for event in events if event.event == "item_succeeded"} == {"rec-a", "rec-b"}

    with repo.connection() as conn:
        runs = conn.execute("SELECT receptor_profile_id, is_current, raw.relative_path FROM docking_runs d JOIN artifacts raw ON raw.artifact_id=d.raw_output_artifact_id ORDER BY receptor_profile_id").fetchall()
        campaign = conn.execute("SELECT workflow_run_id,status,settings_json FROM workflow_runs WHERE workflow_type='docking' ORDER BY started_at DESC LIMIT 1").fetchone()
        stage = conn.execute("SELECT status,summary_json FROM workflow_run_stages WHERE workflow_run_id=? AND stage_name='vina'", (campaign["workflow_run_id"],)).fetchone()
        linked = conn.execute("SELECT COUNT(*) FROM docking_runs WHERE workflow_run_id=?", (campaign["workflow_run_id"],)).fetchone()[0]
    assert [(run["receptor_profile_id"], run["is_current"]) for run in runs] == [("rec-a", 1), ("rec-b", 1)]
    assert all(f"artifacts/docking/{run['receptor_profile_id']}/" in run["relative_path"] for run in runs)
    assert campaign["status"] == "completed" and '"requested_dockings": 2' in campaign["settings_json"]
    assert stage["status"] == "completed" and '"newly_completed": 2' in stage["summary_json"]
    assert linked == 2

    def fake_export_sdf(self, input_pdbqt, output_sdf, log):
        output_sdf.parent.mkdir(parents=True, exist_ok=True)
        output_sdf.write_text("SDF", encoding="ascii")
        log.write_text("exported", encoding="ascii")

    monkeypatch.setattr("moldockpipe.pipeline.PostDockService.export_sdf", fake_export_sdf)
    config = repo.get_settings()
    config["postdock"] = {"mode": "sdf_only", "poses_per_compound": 3, "selected_parents": ["lig1"]}
    (repo.root / "project.yml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    assert runner.run_postdock() == (2, 0)
    for profile in profiles:
        exported = repo.root / "For_PostDocking" / profile["id"] / f"{profile['name']}_prepared.pdb"
        assert exported.is_file()
        assert (exported.parent / "SDF" / "lig1").is_dir()
        assert (exported.parent / "reference_ligand_heavy.mol2").read_text(encoding="utf-8") == standard_mol2

    # Re-running an already exported project refreshes receptor-level package
    # files without re-exporting all ligand poses.
    first_root = repo.root / "For_PostDocking" / profiles[0]["id"]
    (first_root / f"{profiles[0]['name']}_prepared.pdb").unlink()
    (first_root / "reference_ligand_heavy.mol2").unlink()
    assert runner.run_postdock() == (0, 0)
    assert (first_root / f"{profiles[0]['name']}_prepared.pdb").is_file()
    assert (first_root / "reference_ligand_heavy.mol2").is_file()

    calls.clear()
    assert runner.run_vina() == (2, 0)
    assert calls == []

    profiles[1]["exhaustiveness"] = 12
    repo.save_receptor_profiles(profiles)
    assert runner.run_vina({"rec-b"}) == (1, 0)
    assert calls == ["rec-b"]
    with repo.connection() as conn:
        assert conn.execute("SELECT COUNT(*) FROM docking_runs WHERE receptor_profile_id='rec-a'").fetchone()[0] == 1
        assert conn.execute("SELECT COUNT(*) FROM docking_runs WHERE receptor_profile_id='rec-b'").fetchone()[0] == 2
        assert conn.execute("SELECT COUNT(*) FROM docking_runs WHERE receptor_profile_id='rec-b' AND is_current=1").fetchone()[0] == 1

    manifest_outputs = export_manifest_csv(repo)
    leaderboard_outputs = export_leaderboard_csv(repo)
    assert len(manifest_outputs) == len(leaderboard_outputs) == 3
    with (repo.root / "exports" / "manifest.csv").open(encoding="utf-8", newline="") as handle:
        manifest_rows = list(csv.DictReader(handle))
    with (repo.root / "exports" / "leaderboard.csv").open(encoding="utf-8", newline="") as handle:
        leaderboard_rows = list(csv.DictReader(handle))
    assert {(row["receptor_id"], row["id"]) for row in manifest_rows} == {("rec-a", "lig1"), ("rec-b", "lig1")}
    assert {(row["receptor_id"], row["parent_id"]) for row in leaderboard_rows} == {("rec-a", "lig1"), ("rec-b", "lig1")}
