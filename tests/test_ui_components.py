from __future__ import annotations

from types import SimpleNamespace

import pytest

from moldockpipe.project import ProjectRepository

try:
    from moldockpipe.ui.compound_selector import CompoundSelectorDialog
    from moldockpipe.ui.main_window import MainWindow
    from moldockpipe.ui.progress import CheckpointProgress
    from moldockpipe.ui.receptor_manager import ReceptorManagerDialog, ReceptorProfileDialog
    from moldockpipe.ui.receptor_wizard import ReceptorPreparationWizard, ReceptorPreparationWorker
    from moldockpipe.ui.redocking_dialog import RedockingProgressDialog, RedockingResultDialog, RedockingSetupDialog, RedockingWorker
    from moldockpipe.ui.redocking_queue import RedockingQueueController, RedockingQueueDialog
    from moldockpipe.ui.results_dialog import DockingResultsDialog
    from moldockpipe.ui.settings_dialog import SettingsDialog
    from moldockpipe.ui.workers import PipelineWorker
    from moldockpipe.ui.troubleshooter import TroubleshooterDialog
    import moldockpipe.ui.workers as workers
except ImportError as exc:
    pytest.skip(f"Qt UI modules unavailable: {exc}", allow_module_level=True)


def test_ui_components_import_and_expose_expected_contracts() -> None:
    assert SettingsDialog.workflow_reset is not None
    assert CompoundSelectorDialog.selected
    assert CheckpointProgress.STAGES == ("screening", "molscrub", "meeko", "vina", "postdock")
    assert ReceptorManagerDialog is not None
    assert ReceptorProfileDialog is not None
    assert ReceptorPreparationWizard is not None
    assert ReceptorPreparationWorker is not None
    assert RedockingSetupDialog is not None
    assert RedockingProgressDialog is not None
    assert RedockingResultDialog is not None
    assert RedockingWorker is not None
    assert RedockingQueueController is not None
    assert RedockingQueueDialog is not None
    assert DockingResultsDialog is not None
    assert TroubleshooterDialog is not None


def test_troubleshooter_uses_open_project_or_application_root(monkeypatch, tmp_path):
    import moldockpipe.ui.troubleshooter as troubleshooter
    from moldockpipe.services import executables

    roots = []

    class Dialog:
        def __init__(self, root, parent):
            roots.append(root)

        def exec(self):
            pass

    monkeypatch.setattr(troubleshooter, "TroubleshooterDialog", Dialog)
    monkeypatch.setattr(executables, "application_root", lambda: tmp_path / "app")
    MainWindow.show_troubleshooter(SimpleNamespace(repo=None))
    MainWindow.show_troubleshooter(SimpleNamespace(repo=SimpleNamespace(root=tmp_path / "project")))
    assert roots == [tmp_path / "app", tmp_path / "project"]


def test_pipeline_worker_routes_requested_stage(monkeypatch, tmp_path) -> None:
    calls: list[str] = []

    class FakeRunner:
        def __init__(self, repository, progress=None, cancelled=None) -> None:
            self.repository = repository

        def run_screening(self):
            calls.append("screening")
            return (1, 0)

        def run_molscrub(self):
            calls.append("molscrub")
            return (1, 0)

        def run_meeko(self):
            calls.append("meeko")
            return (1, 0)

        def run_vina(self):
            calls.append("vina")
            return (1, 0)

        def run_postdock(self):
            calls.append("postdock")
            return (1, 0)

    monkeypatch.setattr(workers, "PipelineRunner", FakeRunner)
    summaries = []
    worker = PipelineWorker(SimpleNamespace(root=tmp_path, import_project_input=lambda: calls.append("import")), ["meeko"])
    worker.finished.connect(summaries.append)
    worker.run()

    assert calls == ["import", "meeko"]
    assert summaries == [{"meeko": (1, 0)}]


def test_run_all_checks_receptor_before_import_or_chemistry(tmp_path):
    repo = ProjectRepository.create(tmp_path / "project")
    (repo.root / "inputs" / "input.csv").write_text("id,smiles\nlig1,CCO\n", encoding="utf-8")
    worker = PipelineWorker(repo, [])
    failures = []
    worker.failed.connect(failures.append)
    worker.run()
    assert len(failures) == 1 and "Open Receptors" in failures[0]
    with repo.connection() as conn:
        assert conn.execute("SELECT COUNT(*) FROM parent_ligands").fetchone()[0] == 0
        assert conn.execute("SELECT COUNT(*) FROM stage_runs").fetchone()[0] == 0


def test_stale_receptor_manager_checks_before_file_or_database_changes(tmp_path):
    repo = ProjectRepository.create(tmp_path / "project")
    original = {"id": "original", "name": "Original", "enabled": True}
    repo.save_receptor_profiles([original])
    fake = SimpleNamespace(repository=repo, original_profiles={"original": dict(original)})
    newer = {"id": "newer", "name": "Newer", "enabled": True}
    repo.save_receptor_profiles([original, newer])
    with pytest.raises(RuntimeError, match="changed since"):
        ReceptorManagerDialog._save_locked(fake)
    assert repo.get_receptor_profiles() == [original, newer]


def test_bulk_progress_is_throttled_without_losing_failures_or_boundaries(monkeypatch):
    from moldockpipe.models import ProgressEvent

    monkeypatch.setattr(workers, "monotonic", lambda: 10.0)
    worker = PipelineWorker(SimpleNamespace())
    delivered = []
    worker.progress.connect(delivered.append)
    start = ProgressEvent("stage_started", "screening", total=5000)
    worker._report_progress(start)
    for index in range(1, 5001):
        worker._report_progress(ProgressEvent("item_started", "screening", str(index), index, 5000))
        worker._report_progress(ProgressEvent("item_succeeded", "screening", str(index), index, 5000))
    failure = ProgressEvent("item_failed", "screening", "bad", 100, 5000, failed=1)
    worker._report_progress(failure)
    cached_failure = ProgressEvent("item_skipped", "screening", "bad", 101, 5000, failed=1)
    worker._report_progress(cached_failure)
    finish = ProgressEvent("stage_completed", "screening", total=5000)
    worker._report_progress(finish)
    assert delivered[0] == start
    assert [(e.event, e.index) for e in delivered[1:3]] == [
        ("item_started", 5000), ("item_succeeded", 5000)]
    assert delivered[3:] == [failure, cached_failure, finish]
    monkeypatch.setattr(workers, "monotonic", lambda: 10.2)
    update = ProgressEvent("item_succeeded", "meeko", "next", 1, 5000)
    worker._report_progress(update)
    assert delivered[-1] == update
    planned = ProgressEvent("planning_completed", "meeko", total=20, skipped=4980)
    worker._report_progress(planned)
    assert delivered[-1] == planned


def test_dashboard_queries_are_bounded_during_bulk_progress(monkeypatch):
    from moldockpipe.models import ProgressEvent
    import moldockpipe.ui.main_window as main_window

    class Sink:
        def setText(self, value): self.value = value
        def text(self): return getattr(self, "value", "")
        def setValue(self, value): pass
        def setFormat(self, value): pass
        def showMessage(self, value): pass

    refreshes, checkpoints = [], []
    fake = SimpleNamespace(
        stage_progress=Sink(), current_stage_label=Sink(), current_item_label=Sink(),
        current_task_label=Sink(), dashboard_status=Sink(), stat_cards={"running": Sink()},
        stage_started_at=0.0, _refresh_live_dashboard=lambda: refreshes.append(1),
        statusBar=lambda: Sink(), refresh_tables=lambda: None,
        _refresh_checkpoint_state=lambda: checkpoints.append(1),
        overall_progress=SimpleNamespace(complete_stage=lambda stage: None),
        _write_log=lambda message: None, _event_summary=lambda event: "finished",
    )
    monkeypatch.setattr(main_window.time, "monotonic", lambda: 10.0)
    for index in range(1, 5001):
        MainWindow._on_pipeline_progress(fake, ProgressEvent("item_succeeded", "screening", str(index), index, 5000))
    assert len(refreshes) == 1
    monkeypatch.setattr(main_window.time, "monotonic", lambda: 11.1)
    MainWindow._on_pipeline_progress(fake, ProgressEvent("item_succeeded", "screening", "next", 1, 5000))
    assert len(refreshes) == 2
    MainWindow._on_pipeline_progress(fake, ProgressEvent("stage_completed", "screening", total=5000))
    assert len(refreshes) == 3
    assert checkpoints == [1]


def test_docking_completion_distinguishes_new_work_from_reuse():
    from moldockpipe.models import ProgressEvent

    event = ProgressEvent("stage_completed", "vina", succeeded=5000, skipped=4998, failed=1)
    assert MainWindow._event_summary(event) == "Docking completed: newly docked=2, reused=4998, failed=1"


def test_finishing_thread_still_blocks_a_new_pipeline():
    fake = SimpleNamespace(stage_running=False, pipeline_thread=SimpleNamespace(isRunning=lambda: True))
    assert MainWindow._pipeline_active(fake)
    fake.repo = object()
    fake._pipeline_active = lambda: MainWindow._pipeline_active(fake)
    # Must return before constructing a new QThread or changing the dashboard.
    MainWindow._start_pipeline(fake, ["molscrub"], "MolScrub")


@pytest.mark.parametrize("action", ["create_project", "open_project", "import_csv"])
def test_project_changes_are_blocked_before_opening_dialogs(action):
    fake = SimpleNamespace(_block_project_change=lambda: True)
    getattr(MainWindow, action)(fake)


def test_close_waits_for_pipeline_thread_cleanup():
    stopped, ignored, closed = [], [], []
    fake = SimpleNamespace(_pipeline_active=lambda: True, stop_pipeline=lambda: stopped.append(1),
                           close=lambda: closed.append(1), pipeline_thread=object(), pipeline_worker=object())
    MainWindow.closeEvent(fake, SimpleNamespace(ignore=lambda: ignored.append(1)))
    assert ignored == stopped == [1]
    assert closed == []
    assert fake._close_after_pipeline
    MainWindow._pipeline_thread_finished(fake)
    assert fake.pipeline_thread is None and fake.pipeline_worker is None
    assert closed == [1]


def test_dashboard_refresh_defines_profile_metadata(tmp_path) -> None:
    class Sink:
        def setText(self, value): self.text = value
        def setEnabled(self, value): self.enabled = value
        def setStyleSheet(self, value): pass
        def setToolTip(self, value): pass

    fake = SimpleNamespace(
        repo=ProjectRepository.create(tmp_path / "project"), stage_running=False,
        stat_cards={key: Sink() for key in ("ligands", "passed", "failed", "running")},
        stat_details={key: Sink() for key in ("passed", "failed")},
        project_summary=Sink(), docked_export_action=Sink(),
        stage_buttons={key: Sink() for key in ("screening", "molscrub", "meeko", "vina", "postdock")},
        overall_progress=SimpleNamespace(states={}),
    )

    MainWindow._update_dashboard(fake, 0, 0, 0, 0, 0)

    assert fake.docked_export_action.enabled is False

    # Legacy failure notes must remain visible after reopening a project.
    csv_path = tmp_path / "input.csv"
    csv_path.write_text("id,smiles\nlig1,CCO\nlig2,CC\n", encoding="utf-8")
    fake.repo.import_ligands_csv(csv_path)
    with fake.repo.connection() as conn:
        conn.execute("UPDATE parent_ligands SET parse_reason='MolScrub failed: ETKDG failed'")
    MainWindow._update_dashboard(fake, 2, 0, 0, 0, 0)
    assert fake.stat_cards["failed"].text == "2"
    assert "0 generated; 2 failed" in fake.stage_buttons["molscrub"].text
    assert "Failed 3D Structure Generation: 2" in fake.stat_details["failed"].text

    class Progress:
        def reset(self): self.states = {}
        def update(self): pass
        def complete_stage(self, stage): self.states[stage] = "complete"
        def start_stage(self, stage): self.states[stage] = "active"

    fake.overall_progress = Progress()
    fake._update_dashboard = lambda *values: MainWindow._update_dashboard(fake, *values)
    MainWindow._refresh_checkpoint_state(fake)
    assert fake.overall_progress.states["molscrub"] == "warning"
    # An incomplete batch must not become complete after the first result.
    with fake.repo.connection() as conn:
        conn.execute("UPDATE parent_ligands SET parse_reason=NULL WHERE parent_id='lig2'")
    MainWindow._refresh_checkpoint_state(fake)
    assert fake.overall_progress.states["molscrub"] == "active"


def test_auto_import_skips_empty_starter_csv(tmp_path) -> None:
    repo = ProjectRepository.create(tmp_path / "project")
    fake = SimpleNamespace(
        repo=repo,
        _project_input_csv=lambda: repo.root / "inputs" / "input.csv",
    )

    MainWindow._auto_import_input(fake)

    with repo.connection() as conn:
        assert conn.execute("SELECT COUNT(*) FROM parent_ligands").fetchone()[0] == 0


def test_input_csv_has_ligand_rows_excludes_header_only_file(tmp_path) -> None:
    path = tmp_path / "input.csv"
    path.write_text("id,smiles,notes,params_json\n", encoding="utf-8")

    assert MainWindow._input_csv_has_ligand_rows(path) is False

    path.write_text("id,smiles,notes,params_json\nlig1,CCO,,\n", encoding="utf-8")
    assert MainWindow._input_csv_has_ligand_rows(path) is True
