import pytest

from moldockpipe.services.molscrub import MolScrubService
from moldockpipe.pipeline import PipelineRunner
from moldockpipe.project import ProjectRepository


def test_non_enumerating_generation_preserves_a_single_state() -> None:
    pytest.importorskip("rdkit")
    pytest.importorskip("molscrub")
    states, truncated = MolScrubService().generate("CCO", enumerate_states=False, max_states=4)
    assert len(states) == 1
    assert not truncated
    assert states[0].state_isomeric_smiles
    assert states[0].structure_hash


@pytest.mark.parametrize("failure", ["exception", "empty", "persistence"])
def test_molscrub_failure_is_recorded_and_other_ligands_continue(tmp_path, monkeypatch, failure):
    pytest.importorskip("rdkit")
    pytest.importorskip("molscrub")
    repo = ProjectRepository.create(tmp_path / "project")
    csv_path = tmp_path / "input.csv"
    csv_path.write_text("id,smiles\nbad,CC\ngood,CCO\n", encoding="utf-8")
    repo.import_ligands_csv(csv_path)
    original = MolScrubService.generate

    def generate(self, smiles, **kwargs):
        if smiles == "CC":
            if failure == "exception":
                raise RuntimeError("ETKDG conformer generation failed")
            if failure == "empty":
                return [], False
        return original(self, smiles, enumerate_states=False)

    monkeypatch.setattr(MolScrubService, "generate", generate)
    events = []
    runner = PipelineRunner(repo, progress=events.append)
    persist = runner._persist_states

    def persist_then_fail(parent_id, *args):
        persist(parent_id, *args)
        if parent_id == "bad":
            raise RuntimeError("SDF persistence failed")

    if failure == "persistence":
        monkeypatch.setattr(runner, "_persist_states", persist_then_fail)
    assert runner.run_molscrub() == (1, 1)
    assert [row["parent_id"] for row in repo.molscrub_failures()] == ["bad"]
    with repo.connection() as conn:
        assert conn.execute("SELECT parent_id FROM molecular_states WHERE active=1").fetchall()[0][0] == "good"
        assert conn.execute("SELECT status FROM stage_runs WHERE entity_id='bad'").fetchone()[0] == "failed"
    assert events[-1].failed == 1
    assert any(event.event == "item_failed" and event.item_id == "bad" for event in events)

    def unexpected_generation(*args, **kwargs):
        pytest.fail("Unchanged successes and failures must both be reused")

    monkeypatch.setattr(MolScrubService, "generate", unexpected_generation)
    events.clear()
    assert runner.run_molscrub() == (1, 1)
    assert any(event.event == "item_skipped" and event.item_id == "bad" and event.failed == 1 for event in events)
    with repo.connection() as conn:
        assert conn.execute("SELECT COUNT(*) FROM stage_runs WHERE stage_name='molscrub'").fetchone()[0] == 2

    # A successful retry supersedes the failure and clears the legacy note.
    monkeypatch.setattr(MolScrubService, "generate", lambda self, smiles, **kwargs: original(self, smiles, enumerate_states=False))
    monkeypatch.setattr(runner, "_persist_states", persist)
    assert runner.run_molscrub(retry_failed=True) == (2, 0)
    assert repo.molscrub_failures() == []
    with repo.connection() as conn:
        assert conn.execute("SELECT parse_reason FROM parent_ligands WHERE parent_id='bad'").fetchone()[0] is None


def test_isolated_generation_preserves_3d_coordinates():
    pytest.importorskip("rdkit")
    pytest.importorskip("molscrub")
    states, truncated = MolScrubService().generate("CCO", timeout_seconds=15)
    assert states
    assert not truncated
    assert states[0].molecule.GetNumConformers() > 0
    assert states[0].molecule.GetConformer().Is3D()


def test_generation_timeout_stops_worker():
    pytest.importorskip("rdkit")
    with pytest.raises(TimeoutError, match="exceeded.*skipped"):
        MolScrubService().generate("CCO", timeout_seconds=0.001)


def test_disconnected_input_is_rejected_without_starting_worker(monkeypatch):
    pytest.importorskip("rdkit")
    monkeypatch.setattr("moldockpipe.services.molscrub.subprocess.run",
                        lambda *args, **kwargs: pytest.fail("Invalid fragments should fail before worker startup"))
    with pytest.raises(ValueError, match="Disconnected structure"):
        MolScrubService().generate("CCO.[Na+]")
