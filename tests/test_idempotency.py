import sqlite3
import json
from types import SimpleNamespace

import pytest
import yaml

from moldockpipe.pipeline import PipelineRunner
from moldockpipe.project import ProjectRepository


@pytest.fixture
def docking_cache_project(tmp_path, monkeypatch):
    from moldockpipe.project import DEFAULT_VINA_PROFILE

    repo = ProjectRepository.create(tmp_path / "project")
    source = repo.root / "inputs" / "input.csv"
    source.write_text("id,smiles\nlig1,CCO\n", encoding="utf-8")
    repo.import_project_input()
    sdf = repo.root / "state.sdf"
    sdf.write_text("SDF", encoding="ascii")
    sdf_id = repo.add_artifact(sdf, "state_sdf", "molscrub")
    ligand = repo.root / "artifacts" / "pdbqt" / "lig1" / "state1" / "ligand.pdbqt"
    ligand.parent.mkdir(parents=True)
    ligand.write_text("LIGAND", encoding="ascii")
    repo.add_artifact(ligand, "ligand_pdbqt", "meeko")
    with repo.connection() as conn:
        conn.execute("""INSERT INTO molecular_states
            (state_id,parent_id,state_smiles,state_isomeric_smiles,formal_charge,
             state_structure_hash,generation_fingerprint,status,enumeration_truncated,created_at)
            VALUES ('state1','lig1','CCO','CCO',0,'structure','generation','pdbqt_ready',0,'now')""")
        conn.execute("INSERT INTO conformers VALUES ('conf1','state1',1,?,'pdbqt_ready','now')", (sdf_id,))
    receptor = repo.root / "receptor.pdbqt"
    receptor.write_text("RECEPTOR", encoding="ascii")
    profile = {**DEFAULT_VINA_PROFILE, "id": "rec1", "name": "Receptor", "enabled": True,
               "receptor": "receptor.pdbqt"}
    repo.save_receptor_profiles([profile])
    executable = repo.root / "vina.exe"
    executable.write_bytes(b"vina")
    monkeypatch.setattr("moldockpipe.pipeline.find_vina_executable", lambda root: executable)
    calls = []

    class Backend:
        def run(self, *, output, log, **kwargs):
            calls.append(1)
            output.parent.mkdir(parents=True)
            output.write_text("POSE", encoding="ascii")
            log.write_text("LOG", encoding="ascii")
            return SimpleNamespace(return_code=0, poses=[(1,-7,0,0)], command=())

    monkeypatch.setattr("moldockpipe.pipeline.VinaDockingBackend", Backend)
    return SimpleNamespace(repo=repo, runner=PipelineRunner(repo), calls=calls,
                           ligand=ligand, receptor=receptor, executable=executable, profile=profile)


@pytest.mark.parametrize("change", ["noncurrent", "missing_log", "metadata", "numeric_types"])
def test_completed_docking_is_reused_without_scientific_input_changes(docking_cache_project, change):
    p = docking_cache_project
    assert p.runner.run_vina() == (1, 0)
    with p.repo.connection() as conn:
        original = conn.execute("SELECT run_id FROM docking_runs").fetchone()[0]
        if change == "noncurrent":
            conn.execute("UPDATE docking_runs SET is_current=0")
        elif change == "missing_log":
            relative = conn.execute("SELECT a.relative_path FROM artifacts a JOIN docking_runs d ON d.log_artifact_id=a.artifact_id").fetchone()[0]
            (p.repo.root / relative).unlink()
    if change == "metadata":
        p.profile.update(reference_ligand={"name": "updated label"}, description="annotation")
        p.repo.save_receptor_profiles([p.profile])
    elif change == "numeric_types":
        p.profile["size_x"] = float(p.profile["size_x"])
        p.repo.save_receptor_profiles([p.profile])
    events = []
    p.runner.progress = events.append
    assert p.runner.run_vina() == (1, 0)
    assert p.calls == [1]
    assert events[-1].skipped == 1
    with p.repo.connection() as conn:
        assert [tuple(row) for row in conn.execute("SELECT run_id,is_current FROM docking_runs")] == [(original, 1)]


@pytest.mark.parametrize("change", ["ligand", "receptor", "box", "tool", "pose", "pose_metadata"])
def test_changed_scientific_inputs_or_damaged_pose_require_docking(docking_cache_project, change):
    p = docking_cache_project
    assert p.runner.run_vina() == (1, 0)
    if change == "box":
        p.profile["size_x"] += 1
        p.repo.save_receptor_profiles([p.profile])
    elif change == "pose":
        with p.repo.connection() as conn:
            relative = conn.execute("SELECT a.relative_path FROM artifacts a JOIN docking_runs d ON d.raw_output_artifact_id=a.artifact_id").fetchone()[0]
        (p.repo.root / relative).write_text("CORRUPT", encoding="ascii")
    elif change == "pose_metadata":
        with p.repo.connection() as conn:
            conn.execute("DELETE FROM docking_poses")
    else:
        {"ligand": p.ligand, "receptor": p.receptor, "tool": p.executable}[change].write_bytes(b"CHANGED")
    assert p.runner.run_vina() == (1, 0)
    assert len(p.calls) == 2


def test_reverting_box_reuses_historical_completed_docking(docking_cache_project):
    p = docking_cache_project
    assert p.runner.run_vina() == (1, 0)
    original_size = p.profile["size_x"]
    p.profile["size_x"] += 1
    p.repo.save_receptor_profiles([p.profile])
    assert p.runner.run_vina() == (1, 0)
    p.profile["size_x"] = original_size
    p.repo.save_receptor_profiles([p.profile])
    assert p.runner.run_vina() == (1, 0)
    assert len(p.calls) == 2
    with p.repo.connection() as conn:
        assert conn.execute("SELECT is_current FROM docking_runs ORDER BY rowid").fetchall()[0][0] == 1


@pytest.mark.parametrize("legacy_state", ["current", "noncurrent", "changed_input"])
def test_legacy_docking_reuse_requires_valid_current_inputs(docking_cache_project, legacy_state):
    from moldockpipe.fingerprints import file_sha256, fingerprint

    p = docking_cache_project
    assert p.runner.run_vina() == (1, 0)
    with p.repo.connection() as conn:
        command = json.loads(conn.execute("SELECT command_json FROM docking_runs").fetchone()[0])
        command.pop("inputs")
        command["settings"]["reference_ligand"] = {"name": "old label"}
        legacy_fp = fingerprint(settings=command["settings"], inputs={"receptor": file_sha256(p.receptor)},
                                tool_version=file_sha256(p.executable))
        conn.execute("UPDATE docking_runs SET command_json=?,settings_fingerprint=?", (json.dumps(command), legacy_fp))
        if legacy_state == "noncurrent":
            conn.execute("UPDATE docking_runs SET is_current=0")
    if legacy_state == "changed_input":
        p.ligand.write_text("CHANGED", encoding="ascii")
    assert p.runner.run_vina() == (1, 0)
    assert len(p.calls) == (1 if legacy_state == "current" else 2)


def test_failed_attempt_does_not_hide_matching_completed_docking(docking_cache_project):
    p = docking_cache_project
    assert p.runner.run_vina() == (1, 0)
    with p.repo.connection() as conn:
        original = conn.execute("SELECT run_id FROM docking_runs").fetchone()[0]
        conn.execute("UPDATE docking_runs SET is_current=0")
        conn.execute("""INSERT INTO docking_runs
            (run_id,state_id,receptor_hash,settings_fingerprint,status,is_current,receptor_profile_id)
            VALUES ('failed-attempt','state1','old','old','failed',1,'rec1')""")
    assert p.runner.run_vina() == (1, 0)
    assert p.calls == [1]
    with p.repo.connection() as conn:
        assert conn.execute("SELECT run_id FROM docking_runs WHERE is_current=1").fetchone()[0] == original


def test_increasing_generation_timeout_preserves_completed_states(tmp_path, monkeypatch):
    pytest.importorskip("rdkit")
    pytest.importorskip("molscrub")
    repo = ProjectRepository.create(tmp_path / "project")
    source = repo.root / "inputs" / "input.csv"
    source.write_text("id,smiles\nlig1,CCO\n", encoding="utf-8")
    repo.import_project_input()
    runner = PipelineRunner(repo)
    assert runner.run_molscrub() == (1, 0)
    config = repo.get_settings()
    config["molscrub"]["timeout_seconds"] = 120
    (repo.root / "project.yml").write_text(yaml.safe_dump(config), encoding="utf-8")

    def unexpected_generation(*args, **kwargs):
        pytest.fail("Increasing a time limit must not regenerate successful states")

    monkeypatch.setattr("moldockpipe.pipeline.MolScrubService.generate", unexpected_generation)
    assert runner.run_molscrub() == (1, 0)


@pytest.mark.parametrize("same_coordinates", [True, False])
def test_regenerated_ligand_reuse_depends_on_exact_input(docking_cache_project, monkeypatch, same_coordinates):
    pytest.importorskip("rdkit")
    pytest.importorskip("molscrub")
    pytest.importorskip("meeko")
    p = docking_cache_project
    from moldockpipe.services.molscrub import MolScrubService
    from dataclasses import replace
    from rdkit import Chem

    generate = MolScrubService.generate
    saved = []

    def regenerate(self, *args, **kwargs):
        if not saved:
            saved.append(generate(self, *args, **kwargs))
            return saved[0]
        states, truncated = saved[0]
        if same_coordinates:
            return states, truncated
        changed = []
        for state in states:
            molecule = Chem.Mol(state.molecule)
            point = molecule.GetConformer().GetAtomPosition(0)
            point.x += 0.25
            molecule.GetConformer().SetAtomPosition(0, point)
            changed.append(replace(state, molecule=molecule))
        return changed, truncated

    monkeypatch.setattr(MolScrubService, "generate", regenerate)
    p.repo.clear_workflow_data()
    p.repo.import_project_input()
    p.runner.run_screening()
    assert p.runner.run_molscrub() == (1, 0)
    assert p.runner.run_meeko() == (1, 0)
    assert p.runner.run_vina() == (1, 0)
    with p.repo.connection() as conn:
        relative = conn.execute("SELECT a.relative_path FROM conformers c JOIN artifacts a ON a.artifact_id=c.sdf_artifact_id").fetchone()[0]
    (p.repo.root / relative).write_text("CORRUPT", encoding="ascii")
    assert p.runner.run_molscrub() == (1, 0)
    assert p.runner.run_meeko() == (1, 0)
    assert p.runner.run_vina() == (1, 0)
    assert len(p.calls) == (1 if same_coordinates else 2)


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


@pytest.mark.parametrize("status", ["running", "interrupted", "missing"])
def test_incomplete_state_generation_is_not_reused(tmp_path, monkeypatch, status):
    pytest.importorskip("rdkit")
    pytest.importorskip("molscrub")
    from moldockpipe.services.molscrub import MolScrubService

    repo = ProjectRepository.create(tmp_path / "project")
    source = repo.root / "inputs" / "input.csv"
    source.write_text("id,smiles\nlig1,CCO\n", encoding="utf-8")
    repo.import_project_input()
    runner = PipelineRunner(repo)
    assert runner.run_molscrub() == (1, 0)
    # A process can die after saving valid states but before finishing the
    # parent's complete state set. File validity alone cannot prove completion.
    with repo.connection() as conn:
        if status == "missing":
            conn.execute("DELETE FROM stage_runs WHERE stage_name='molscrub'")
        else:
            conn.execute("UPDATE stage_runs SET status=? WHERE stage_name='molscrub'", (status,))
    original = MolScrubService.generate
    calls = []

    def generate(self, *args, **kwargs):
        calls.append(1)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(MolScrubService, "generate", generate)
    assert runner.run_molscrub() == (1, 0)
    assert calls == [1]
    assert runner.run_molscrub() == (1, 0)
    assert calls == [1]
