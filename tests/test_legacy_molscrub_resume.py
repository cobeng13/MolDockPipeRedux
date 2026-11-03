"""Historical structures must not silently become new campaign inputs."""
import json
from types import SimpleNamespace

import pytest

from moldockpipe.fingerprints import fingerprint, file_sha256
from moldockpipe.pipeline import PipelineRunner
from moldockpipe.project import DEFAULT_VINA_PROFILE, ProjectRepository


@pytest.fixture
def legacy_project(tmp_path, monkeypatch):
    repo = ProjectRepository.create(tmp_path / "project")
    (repo.root / "inputs" / "input.csv").write_text("id,smiles\nlig1,CCO\n", encoding="utf-8")
    repo.import_project_input()
    settings = {k: v for k, v in repo.get_settings()["molscrub"].items() if k != "workers"}
    state_fp = fingerprint(settings=settings, inputs={"state": "structure"}, tool_version="unknown")
    parent_fp = fingerprint(settings=settings, inputs={"smiles": "CCO"}, tool_version="unknown")
    sdf = repo.root / "saved.sdf"
    sdf.write_text("original saved coordinates", encoding="ascii")
    artifact = repo.add_artifact(sdf, "state_sdf", "molscrub")
    ligand = repo.root / "artifacts" / "pdbqt" / "lig1" / "state1" / "ligand.pdbqt"
    ligand.parent.mkdir(parents=True)
    ligand.write_text("original prepared ligand", encoding="ascii")
    repo.add_artifact(ligand, "ligand_pdbqt", "meeko")
    with repo.connection() as conn:
        conn.execute("""INSERT INTO molecular_states
            (state_id,parent_id,state_smiles,state_isomeric_smiles,formal_charge,
             state_structure_hash,generation_fingerprint,status,created_at)
            VALUES ('state1','lig1','CCO','CCO',0,'structure',?,'pdbqt_ready','now')""", (state_fp,))
        conn.execute("INSERT INTO conformers VALUES ('conf1','state1',1,?,'pdbqt_ready','now')", (artifact,))
        conn.execute("INSERT INTO stage_runs VALUES ('old','molscrub','parent','lig1',?,'completed','now','now',NULL)", (parent_fp,))
    receptor = repo.root / "receptor.pdbqt"
    receptor.write_text("receptor", encoding="ascii")
    repo.save_receptor_profiles([{**DEFAULT_VINA_PROFILE, "receptor": "receptor.pdbqt"}])
    executable = repo.root / "vina.exe"
    executable.write_bytes(b"test backend")
    monkeypatch.setattr("moldockpipe.pipeline.find_vina_executable", lambda root: executable)

    class Backend:
        def run(self, *, output, log, **kwargs):
            output.parent.mkdir(parents=True)
            output.write_text("saved docking pose", encoding="ascii")
            log.write_text("log", encoding="ascii")
            return SimpleNamespace(return_code=0, poses=[(1, -7, 0, 0)], command=())

    monkeypatch.setattr("moldockpipe.pipeline.VinaDockingBackend", Backend)
    assert PipelineRunner(repo).run_vina() == (1, 0)
    monkeypatch.setattr("moldockpipe.pipeline.MolScrubService.version", staticmethod(lambda **kwargs: "0.3.1"))
    return SimpleNamespace(repo=repo, sdf=sdf, ligand=ligand, settings=settings, state_fp=state_fp)


def ledger(repo):
    with repo.connection() as conn:
        return {table: [tuple(row) for row in conn.execute(f"SELECT * FROM {table}")]
                for table in ("stage_runs", "molecular_states", "conformers", "docking_runs", "artifacts")}


@pytest.mark.parametrize("history", [True, False])
def test_version_transition_stops_before_any_scientific_writes(legacy_project, history):
    p = legacy_project
    if not history:
        with p.repo.connection() as conn:
            conn.execute("DELETE FROM stage_runs")
    before = ledger(p.repo)
    digest = file_sha256(p.sdf)
    with pytest.raises(RuntimeError, match="Retain verified legacy results"):
        PipelineRunner(ProjectRepository(p.repo.root)).run_molscrub()
    assert ledger(p.repo) == before
    assert file_sha256(p.sdf) == digest


@pytest.mark.parametrize("history", [True, False])
def test_approved_legacy_resume_preserves_full_downstream_campaign(legacy_project, monkeypatch, history):
    p = legacy_project
    if not history:
        with p.repo.connection() as conn:
            conn.execute("DELETE FROM stage_runs")
    p.repo.update_settings({"resume": {"retain_legacy_molscrub": True}})
    before = ledger(p.repo)
    hashes = [file_sha256(path) for path in (p.sdf, p.ligand)]

    def unexpected(*args, **kwargs):
        pytest.fail("Retained campaign inputs must not dispatch chemistry or docking")

    monkeypatch.setattr("moldockpipe.pipeline.MolScrubService.generate", unexpected)
    monkeypatch.setattr("moldockpipe.pipeline.MeekoService.prepare_ligand", unexpected)
    monkeypatch.setattr("moldockpipe.pipeline.VinaDockingBackend.run", unexpected)
    for _ in range(2):
        fresh = PipelineRunner(ProjectRepository(p.repo.root))
        assert fresh.run_molscrub() == (1, 0)
        assert fresh.run_meeko() == (0, 0)
        assert fresh.run_vina() == (1, 0)
    after = ledger(p.repo)
    # A docking campaign may update its bookkeeping, but keeps scientific identity.
    assert after["docking_runs"] == before["docking_runs"]
    assert after["stage_runs"] == before["stage_runs"]
    assert after["molecular_states"] == before["molecular_states"]
    assert after["conformers"] == before["conformers"]
    assert [file_sha256(path) for path in (p.sdf, p.ligand)] == hashes
    with p.repo.connection() as conn:
        event = conn.execute("SELECT data_json FROM provenance_events WHERE event_type='legacy_molscrub_results_retained'").fetchone()
    assert json.loads(event[0])["original_version"] == "unknown"
    assert json.loads(event[0])["runtime_version"] == "0.3.1"


@pytest.mark.parametrize("change", ["missing", "corrupt", "ph", "max_states", "known_version", "interrupted", "partial", "unprepared", "inactive_artifact"])
def test_legacy_choice_does_not_bypass_input_or_completion_checks(legacy_project, change):
    p = legacy_project
    settings = dict(p.settings)
    if change == "missing":
        p.sdf.unlink()
    elif change == "corrupt":
        p.sdf.write_text("changed", encoding="ascii")
    elif change in {"ph", "max_states"}:
        settings[change] += 1
    else:
        with p.repo.connection() as conn:
            if change == "known_version":
                fp = fingerprint(settings=settings, inputs={"state": "structure"}, tool_version="0.2.2")
                conn.execute("UPDATE molecular_states SET generation_fingerprint=?", (fp,))
            elif change == "interrupted":
                conn.execute("UPDATE stage_runs SET status='interrupted'")
            else:
                conn.execute("DELETE FROM stage_runs")
                if change == "unprepared":
                    conn.execute("UPDATE molecular_states SET status='sdf_ready'")
                    conn.execute("UPDATE conformers SET status='ready'")
                else:
                    conn.execute("UPDATE conformers SET status='pdbqt_failed'" if change == "partial" else
                                 "UPDATE artifacts SET active=0 WHERE created_stage='molscrub'")
    runner = PipelineRunner(p.repo)
    with p.repo.connection() as conn:
        rows = conn.execute("SELECT parent_id,source_smiles FROM parent_ligands").fetchall()
    assert runner._molscrub_resume_plan(rows, settings, "0.3.1")[0][2] == "generate"


def test_current_version_without_old_history_reuses_only_complete_states(legacy_project, monkeypatch):
    p = legacy_project
    with p.repo.connection() as conn:
        conn.execute("DELETE FROM stage_runs")
    monkeypatch.setattr("moldockpipe.pipeline.MolScrubService.version", staticmethod(lambda **kwargs: "unknown"))
    assert PipelineRunner(ProjectRepository(p.repo.root)).run_molscrub() == (1, 0)


def test_legacy_failure_requires_choice_and_explicit_retry_remains_available(legacy_project):
    p = legacy_project
    with p.repo.connection() as conn:
        conn.execute("UPDATE stage_runs SET status='failed',reason='old failure'")
        conn.execute("UPDATE molecular_states SET active=0")
    with pytest.raises(RuntimeError, match="Retain verified legacy results"):
        PipelineRunner(p.repo).run_molscrub()
    p.repo.update_settings({"resume": {"retain_legacy_molscrub": True}, "molscrub": {"timeout_seconds": 15}})
    assert PipelineRunner(p.repo).run_molscrub() == (0, 1)
    with p.repo.connection() as conn:
        rows = conn.execute("SELECT parent_id,source_smiles FROM parent_ligands").fetchall()
    runner = PipelineRunner(p.repo)
    assert runner._molscrub_resume_plan(rows, p.settings, "0.3.1", retry_failed=True)[0][2] == "generate"
    assert runner._molscrub_resume_plan(rows, {**p.settings, "timeout_seconds": 60}, "0.3.1")[0][2] == "generate"


def test_legacy_mismatch_is_found_before_an_earlier_pending_job(legacy_project):
    p = legacy_project
    with p.repo.connection() as conn:
        conn.execute("INSERT INTO parent_ligands(parent_id,source_smiles,parse_status,created_at) VALUES ('aaa','CC','pending','now')")
    before = ledger(p.repo)
    with pytest.raises(RuntimeError, match="Retain verified legacy results"):
        PipelineRunner(p.repo).run_molscrub()
    assert ledger(p.repo) == before


@pytest.mark.parametrize("stored,current", [(None, 15), (15, None), (15.0, 15), (15, 15.0), (15.0, 15.0), (15, 60)])
def test_saving_default_timeout_does_not_retry_unchanged_failures(legacy_project, stored, current):
    p = legacy_project
    old_settings = dict(p.settings)
    new_settings = dict(p.settings)
    if stored is not None:
        old_settings["timeout_seconds"] = stored
    if current is not None:
        new_settings["timeout_seconds"] = current
    old_fp = fingerprint(settings=old_settings, inputs={"smiles": "CCO"}, tool_version="0.3.1")
    with p.repo.connection() as conn:
        conn.execute("UPDATE stage_runs SET status='failed',fingerprint=?", (old_fp,))
        rows = conn.execute("SELECT parent_id,source_smiles FROM parent_ligands").fetchall()
    decision = PipelineRunner(p.repo)._molscrub_resume_plan(rows, new_settings, "0.3.1")[0][2]
    assert decision == ("generate" if current == 60 else "failure")
