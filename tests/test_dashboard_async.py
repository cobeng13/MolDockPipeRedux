"""Responsiveness and lifecycle regressions for large-campaign UI work."""
from threading import Event, get_ident
from types import SimpleNamespace

import pytest

pytest.importorskip("PyQt6")
from PyQt6.QtCore import QTimer, QCoreApplication, QEvent
from PyQt6.QtTest import QTest
from test_results_viewer import qt_app, wait_for
from moldockpipe.project import ProjectRepository
from moldockpipe.ui.main_window import MainWindow


def test_slow_dashboard_keeps_gui_alive_and_coalesces_requests(qt_app, tmp_path, monkeypatch):
    import moldockpipe.ui.main_window as module
    from moldockpipe.dashboard import collect_dashboard
    repo = ProjectRepository.create(tmp_path / "project")
    entered, release = Event(), Event()
    calls, ticks = [], []
    stopped = Event()
    gui_thread = get_ident()

    def slow(repo, running, cancelled):
        calls.append(get_ident())
        entered.set()
        assert release.wait(5)
        return collect_dashboard(repo, running, cancelled)

    monkeypatch.setattr(module, "collect_dashboard", slow)
    window = MainWindow()
    window.repo = repo
    window.stage_running = True
    window.pipeline_worker = SimpleNamespace(stop=stopped.set)
    timer = QTimer()
    timer.timeout.connect(lambda: ticks.append(1))
    timer.start(10)
    try:
        window._request_dashboard()
        wait_for(qt_app, entered.is_set)
        for _ in range(20000):
            window._request_dashboard()
        QTest.qWait(100)
        assert len(ticks) >= 3
        assert calls == [calls[0]] and calls[0] != gui_thread
        window.stop_pipeline()
        assert stopped.is_set()
        release.set()
        wait_for(qt_app, lambda: window._dashboard_worker is None)
        assert window._dashboard_timer.remainingTime() > 500
        assert window._dashboard_values == (0, 0, 0, 0, 0)
    finally:
        release.set()
        wait_for(qt_app, lambda: window._dashboard_worker is None)
        timer.stop()
        window.stage_running = False
        window.pipeline_worker = None
        window.close()
        window.deleteLater()
        QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)


def test_stale_dashboard_cannot_overwrite_new_project_or_running_state():
    old, new = object(), object()
    applied = []
    fake = SimpleNamespace(repo=new, stage_running=False, _closing=False,
        _dashboard_pending=False, _apply_dashboard=lambda *args: applied.append(args))
    MainWindow._dashboard_ready(fake, old, False, {})
    MainWindow._dashboard_ready(fake, new, True, {})
    assert not applied and fake._dashboard_pending


def test_close_waits_for_snapshot_without_blocking_gui(qt_app, tmp_path, monkeypatch):
    import moldockpipe.ui.main_window as module
    entered, release = Event(), Event()

    def snapshot(repo, running, cancelled):
        entered.set()
        assert release.wait(5)
        return {}

    monkeypatch.setattr(module, "collect_dashboard", snapshot)
    window = MainWindow()
    window.repo = ProjectRepository.create(tmp_path / "project")
    window.show()
    try:
        window._request_dashboard()
        wait_for(qt_app, entered.is_set)
        window.close()
        assert window._closing and window.isVisible()
        assert window._dashboard_worker.isInterruptionRequested()
        release.set()
        wait_for(qt_app, lambda: window._dashboard_worker is None)
        assert not window.isVisible()
        assert not window._dashboard_timer.isActive()
    finally:
        release.set()
        wait_for(qt_app, lambda: window._dashboard_worker is None)
        window.close()
        window.deleteLater()
        QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)


def test_project_io_runs_off_gui_thread(qt_app, tmp_path):
    window = MainWindow()
    results = []
    gui_thread = get_ident()
    window._start_project_task(lambda: get_ident(), results.append, "Loading")
    wait_for(qt_app, lambda: window._project_task is None)
    assert results and results[0] != gui_thread
    window.close()
    window.deleteLater()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)


def test_performance_report_contains_counts_and_resource_sample(tmp_path):
    import json
    from moldockpipe.performance import save_performance_report
    repo = ProjectRepository.create(tmp_path / "project")
    output = tmp_path / "performance.json"
    save_performance_report(repo, output, {"heartbeat_max_gap_seconds": 0.25})
    report = json.loads(output.read_text())
    assert report['record_counts']['parent_ligands'] == 0
    assert 'idx_stage_runs_entity' in report['indexes']['stage_runs']
    assert report['gui']['heartbeat_max_gap_seconds'] == 0.25
    if __import__('os').name == 'nt':
        assert report['available_memory_bytes'] > 0
        current = next(row for row in report['processes'] if row['pid'] == __import__('os').getpid())
        assert current['private_bytes'] > 0
