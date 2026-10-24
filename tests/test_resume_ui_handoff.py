"""UI resume regressions independent of an interactive Qt desktop."""
from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("PyQt6")
from moldockpipe.models import ProgressEvent
from moldockpipe.project import DEFAULT_VINA_PROFILE, ProjectRepository
from moldockpipe.ui.main_window import MainWindow
from moldockpipe.ui.workers import PipelineWorker
import moldockpipe.ui.workers as workers


class Sink:
    def setText(self, text): self.text = text
    def setEnabled(self, enabled): self.enabled = enabled
    def setStyleSheet(self, value): pass
    def setToolTip(self, value): pass
    def setRange(self, low, high): self.range = (low, high)
    def setValue(self, value): self.value = value
    def setFormat(self, value): self.format = value


@pytest.mark.parametrize("running", [True, False])
def test_active_dashboard_never_walks_export_directories(tmp_path, monkeypatch, running):
    repo = ProjectRepository.create(tmp_path / "project")
    repo.save_receptor_profiles([{**DEFAULT_VINA_PROFILE, "id": "rec", "enabled": True}])
    with repo.connection() as conn:
        conn.execute("INSERT INTO parent_ligands(parent_id,source_smiles,parse_status,created_at) VALUES ('p','CCO','completed','now')")
        conn.execute("""INSERT INTO molecular_states
            (state_id,parent_id,state_smiles,state_isomeric_smiles,formal_charge,
             state_structure_hash,generation_fingerprint,status,created_at)
            VALUES ('s','p','CCO','CCO',0,'hash','fp','pdbqt_ready','now')""")
        conn.execute("""INSERT INTO docking_runs
            (run_id,state_id,receptor_hash,settings_fingerprint,status,is_current,receptor_profile_id)
            VALUES ('run','s','receptor','fp','completed',1,'rec')""")
    from moldockpipe.dashboard import collect_dashboard
    scans = []
    def scanned(path):
        scans.append(path)
        if running:
            pytest.fail("Active dashboard performed filesystem export scan")
        return False
    monkeypatch.setattr(Path, "is_dir", scanned)
    snapshot = collect_dashboard(repo, running)
    assert snapshot['has_docking_exports']
    assert snapshot['export_status'] == 'Ready'
    assert bool(scans) == (not running)


def test_worker_announces_stage_before_dependency_setup(tmp_path, monkeypatch):
    events, calls = [], []

    class Runner:
        def __init__(self, *args, **kwargs): pass
        def run_molscrub(self):
            assert events[-1].event == "stage_preparing"
            assert events[-1].stage == "molscrub"
            calls.append("load dependencies")
            return 1, 0

    monkeypatch.setattr(workers, "PipelineRunner", Runner)
    worker = PipelineWorker(SimpleNamespace(root=tmp_path, import_project_input=lambda: None), ["molscrub"])
    worker.progress.connect(events.append)
    worker.run()
    assert calls == ["load dependencies"]


def test_stop_at_stage_handoff_does_not_enter_dependency_setup(tmp_path, monkeypatch):
    class Runner:
        def __init__(self, *args, **kwargs): pass
        def run_molscrub(self): pytest.fail("Stopped worker entered next stage")

    monkeypatch.setattr(workers, "PipelineRunner", Runner)
    worker = PipelineWorker(SimpleNamespace(root=tmp_path, import_project_input=lambda: None), ["molscrub"])
    worker.progress.connect(lambda event: worker.stop())
    stopped = []
    worker.stopped.connect(lambda: stopped.append(True))
    worker.run()
    assert stopped == [True]


def test_dependency_failure_clears_busy_progress_and_running_labels(monkeypatch):
    import moldockpipe.ui.main_window as window
    monkeypatch.setattr(window.QMessageBox, "critical", lambda *args: None)
    fake = SimpleNamespace(
        stage_running=True, stage_progress=Sink(),
        current_stage_label=Sink(), current_item_label=Sink(), current_task_label=Sink(),
        dashboard_status=Sink(), run_all_action=Sink(), refresh_action=Sink(), stop_action=Sink(),
        _set_stage_actions=lambda enabled: None, _write_log=lambda message: None,
        _refresh_checkpoint_state=lambda: None,
    )
    MainWindow._pipeline_failed(fake, "Loading MolScrub dependencies exceeded 30s")
    assert fake.stage_running is False
    assert fake.stage_progress.range == (0, 100)
    assert fake.dashboard_status.text == "Pipeline failed"
    assert "exceeded 30s" in fake.current_task_label.text
    assert fake.stop_action.enabled is False


@pytest.mark.parametrize('kind', ['item_skipped', 'item_failed'])
def test_failure_bursts_are_rate_limited_without_losing_final_totals(monkeypatch, kind):
    monkeypatch.setattr(workers, 'monotonic', lambda: 0.0)
    worker = PipelineWorker(None)
    events = []
    worker.progress.connect(events.append)
    for index in range(1, 20001):
        worker._report_progress(ProgressEvent(kind, 'molscrub', index=index, total=20000, failed=index))
    worker._report_progress(ProgressEvent('stage_completed', 'molscrub', total=20000, failed=20000))
    assert len(events) == 3
    assert events[-1].failed == 20000
