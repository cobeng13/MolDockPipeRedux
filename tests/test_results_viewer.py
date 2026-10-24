"""UI control/lifecycle tests use a fake web canvas; real rendering is opt-in."""
import base64
import os
import time

import pytest

from test_structure_viewer import viewer_project, PDB, SDF


@pytest.fixture(scope="session")
def qt_app():
    pytest.importorskip("PyQt6")
    from moldockpipe.__main__ import configure_qt_runtime
    configure_qt_runtime()
    # WebEngine requires one QApplication for the entire process; recreating it
    # after a complete teardown is unsupported by Chromium's global runtime.
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PyQt6.QtCore import QCoreApplication, Qt
    from PyQt6.QtWidgets import QApplication
    QCoreApplication.setAttribute(Qt.ApplicationAttribute.AA_ShareOpenGLContexts)
    app = QApplication.instance() or QApplication(["moldockpipe-viewer-test"])
    yield app


@pytest.fixture
def fake_canvas(qt_app, monkeypatch):
    from PyQt6.QtCore import pyqtSignal
    from PyQt6.QtWidgets import QWidget
    import moldockpipe.ui.molecular_viewer as module

    class Canvas(QWidget):
        loadFinished = pyqtSignal(bool)
        renderProcessTerminated = pyqtSignal()

        def page(self):
            return self

    monkeypatch.setattr(module, "_web_view", lambda parent, bridge: Canvas(parent))


def wait_for(qt_app, condition):
    from PyQt6.QtTest import QTest
    deadline = time.monotonic() + 8
    while not condition() and time.monotonic() < deadline:
        qt_app.processEvents()
        QTest.qWait(10)
    assert condition()


def test_result_best_mode_and_pose_switch(qt_app, fake_canvas, viewer_project):
    from moldockpipe.ui.results_dialog import DockingResultsDialog
    repo, service, calls, receptor, raw = viewer_project
    # Replace the class's new service with the prepared service/spy.
    import moldockpipe.ui.results_dialog as module
    original = module.StructureViewerService
    module.StructureViewerService = lambda _: service
    try:
        dialog = DockingResultsDialog(repo)
    finally:
        module.StructureViewerService = original
    dialog.viewer.bridge.ready()
    wait_for(qt_app, lambda: dialog._worker is None)
    assert dialog.viewer._structure.mode_index == 8
    assert dialog.viewer.pose.count() == 2
    assert not dialog.viewer.export.isEnabled()
    dialog.viewer.bridge.rendered(dialog.viewer._token)
    assert dialog.viewer.export.isEnabled()
    dialog.viewer.pose.setCurrentIndex(0)
    wait_for(qt_app, lambda: dialog._worker is None)
    assert dialog.viewer._structure.mode_index == 3
    assert "Mode 3" in dialog.viewer.title.text()
    dialog.reject()


def test_pending_row_and_tab_switch_clear_old_pose(qt_app, fake_canvas, viewer_project, monkeypatch):
    from moldockpipe.ui.results_dialog import DockingResultsDialog
    import moldockpipe.ui.results_dialog as module
    repo, service, calls, receptor, raw = viewer_project
    csv = repo.root / "pending.csv"
    csv.write_text("id,smiles\npending,CCC\n")
    repo.import_ligands_csv(csv)
    profiles = repo.get_receptor_profiles()
    repo.save_receptor_profiles(profiles + [{**profiles[0], "id": "other", "name": "Other"}])
    monkeypatch.setattr(module, "StructureViewerService", lambda _: service)
    dialog = DockingResultsDialog(repo)
    dialog.viewer.bridge.ready()
    wait_for(qt_app, lambda: dialog._worker is None)
    table = dialog.tabs.currentWidget()
    pending_row = next(row for row in range(table.rowCount()) if table.item(row, 0).text() == "pending")
    table.selectRow(pending_row)
    assert dialog.viewer._structure is None
    assert "pending" in dialog.viewer.title.text()
    assert not dialog.viewer.export.isEnabled()
    dialog.tabs.setCurrentIndex(1)
    assert dialog.viewer._structure is None
    assert dialog.viewer.pose.count() == 0
    dialog.reject()


def test_stale_load_is_ignored_and_closing_cancels_worker(qt_app, fake_canvas, viewer_project, monkeypatch):
    from moldockpipe.ui.results_dialog import DockingResultsDialog
    import moldockpipe.ui.results_dialog as module
    repo, service, calls, receptor, raw = viewer_project
    started = __import__("threading").Event()
    original = service.load

    def delayed(run_id, pose_id, *, cancelled):
        started.set()
        while not cancelled():
            time.sleep(0.01)
        raise InterruptedError("cancelled")

    monkeypatch.setattr(service, "load", delayed)
    monkeypatch.setattr(module, "StructureViewerService", lambda _: service)
    dialog = DockingResultsDialog(repo)
    wait_for(qt_app, started.is_set)
    previous_generation = dialog._generation
    dialog._invalidate("New selection")
    dialog._loaded(previous_generation, object())  # An old result must never reach the canvas.
    assert dialog.viewer._structure is None
    dialog.reject()
    wait_for(qt_app, lambda: dialog._worker is None)
    assert dialog._closing


def test_latest_selection_runs_after_cancelled_conversion(qt_app, fake_canvas, viewer_project, monkeypatch):
    from moldockpipe.ui.results_dialog import DockingResultsDialog
    import moldockpipe.ui.results_dialog as module
    repo, service, calls, receptor, raw = viewer_project
    started = __import__("threading").Event()
    original = service.load

    def delayed(run_id, pose_id, *, cancelled):
        if not started.is_set():
            started.set()
            while not cancelled():
                time.sleep(0.01)
            raise InterruptedError("cancelled")
        return original(run_id, pose_id, cancelled=cancelled)

    monkeypatch.setattr(service, "load", delayed)
    monkeypatch.setattr(module, "StructureViewerService", lambda _: service)
    dialog = DockingResultsDialog(repo)
    wait_for(qt_app, started.is_set)
    dialog.viewer.pose.setCurrentIndex(0)
    wait_for(qt_app, lambda: dialog._worker is None and dialog._pending is None)
    assert dialog.viewer._structure.mode_index == 3
    dialog.reject()


def test_missing_webengine_keeps_results_table(qt_app, viewer_project, monkeypatch):
    import moldockpipe.ui.molecular_viewer as module
    from moldockpipe.ui.results_dialog import DockingResultsDialog
    monkeypatch.setattr(module, "_web_view", lambda *_: (_ for _ in ()).throw(ImportError("missing WebEngine")))
    dialog = DockingResultsDialog(viewer_project[0])
    assert dialog.tabs.currentWidget().rowCount() == 1
    assert not dialog.viewer.available
    assert "README" in dialog.viewer.status.text()
    assert dialog._worker is None
    dialog.reject()


def test_contact_controls_send_visibility_and_clear_on_pose_change(qt_app, fake_canvas):
    import json
    from moldockpipe.services.structure_viewer import ViewerStructure
    from moldockpipe.ui.molecular_viewer import MolecularViewer
    viewer = MolecularViewer()
    messages = []
    viewer.bridge.command.connect(lambda value: messages.append(json.loads(value)))
    viewer.bridge.ready()
    contact = dict(kind="hbond", start=dict(x=0, y=0, z=0), end=dict(x=0, y=0, z=3),
                   distance=3, residue="ALA A:10", protein_serials=[1], detail="H-bond candidate")
    viewer.set_structure(ViewerStructure(PDB, SDF, "Target", "ligand", "state", "run", "pose", 1, -5,
                                         [contact], "Geometric predictions"))
    assert messages[-1]["interactions"] == [contact]
    assert "hbond" in messages[-1]["interaction_types"]
    assert "hydrophobic" not in messages[-1]["interaction_types"]
    assert viewer.interaction_toggles["hbond"].text() == "H-bonds (1)"
    viewer.bridge.rendered(viewer._token)
    viewer.interaction_toggles["hbond"].setChecked(False)
    assert "hbond" not in messages[-1]["interaction_types"]
    viewer.interaction_labels.setChecked(False)
    assert messages[-1]["interaction_labels"] is False
    viewer.show_receptor.setChecked(False)
    assert messages[-1]["show_receptor"] is False
    viewer.clear("Loading next pose")
    assert messages[-1]["kind"] == "clear"
    assert not viewer.interaction_info.text()
    assert not viewer.interaction_labels.isEnabled()
    viewer.close()


def test_png_save_stale_export_and_write_failure(qt_app, fake_canvas, tmp_path, monkeypatch):
    from PyQt6.QtCore import QBuffer, QIODevice
    from PyQt6.QtGui import QImage, QColor
    from moldockpipe.ui.molecular_viewer import MolecularViewer
    from PyQt6.QtWidgets import QMessageBox
    viewer = MolecularViewer()
    image = QImage(32, 24, QImage.Format.Format_RGB32)
    image.fill(QColor("white"))
    buffer = QBuffer()
    buffer.open(QIODevice.OpenModeFlag.WriteOnly)
    image.save(buffer, "PNG")
    data = "data:image/png;base64," + base64.b64encode(bytes(buffer.data())).decode()
    target = tmp_path / "view.png"
    viewer._export_path = str(target)
    viewer._save_png(viewer._token + 1, data)
    assert not target.exists()
    viewer._save_png(viewer._token, data)
    assert QImage(str(target)).size() == image.size()
    warnings = []
    monkeypatch.setattr(QMessageBox, "warning", lambda *args: warnings.append(args))
    viewer._export_path = str(tmp_path / "missing" / "view.png")
    viewer._save_png(viewer._token, data)
    assert warnings and "failed" in viewer.status.text()
    viewer._export_path = str(target)
    viewer._save_png(viewer._token, "data:image/png;base64,bm90IGEgcG5n")
    assert len(warnings) == 2
    assert QImage(str(target)).size() == image.size()
    viewer.close()


@pytest.mark.skipif(os.environ.get("MOLDOCKPIPE_WEBENGINE_TEST") != "1", reason="Opt-in native WebEngine/WebGL integration")
def test_native_viewer_destroys_page_before_profile(qt_app):
    from PyQt6 import sip
    from PyQt6.QtCore import qInstallMessageHandler
    from PyQt6.QtWidgets import QWidget
    from moldockpipe.ui.molecular_viewer import MolecularViewer
    messages = []
    previous = qInstallMessageHandler(lambda kind, context, message: messages.append(message))
    try:
        for _ in range(3):
            parent = QWidget()
            viewer = MolecularViewer(parent)
            parent.show()
            wait_for(qt_app, lambda: viewer._ready or not viewer.available)
            assert viewer._ready, viewer.status.text()
            profile, _, page, _ = viewer.web._viewer_objects
            destroyed = []
            page.destroyed.connect(lambda: destroyed.append("page"))
            profile.destroyed.connect(lambda: destroyed.append("profile"))
            # Parent deletion also covers application/window shutdown without a
            # viewer closeEvent, rather than relying on queued cleanup callbacks.
            sip.delete(parent)
            assert sip.isdeleted(page) and sip.isdeleted(profile)
            assert destroyed == ["page", "profile"]
            qt_app.processEvents()
    finally:
        qInstallMessageHandler(previous)
    assert not any("WebEnginePage still not deleted" in message for message in messages), messages


@pytest.mark.skipif(os.environ.get("MOLDOCKPIPE_WEBENGINE_TEST") != "1", reason="Opt-in native WebEngine/WebGL integration")
def test_native_offline_render_and_double_size_png(qt_app, tmp_path):
    from PyQt6.QtGui import QImage
    from moldockpipe.ui.molecular_viewer import MolecularViewer
    from moldockpipe.services.structure_viewer import ViewerStructure
    # Connected atoms make a visible ligand; large text exercises the bridge
    # beyond QWebEngineView.setHtml's 2MB limit without an oversized atom graph.
    ligand = """ethanol
  MolDockPipe

  3  2  0  0  0  0            999 V2000
    1.0000    2.0000    3.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5000    2.0000    3.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2000    3.2000    3.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$
"""
    viewer = MolecularViewer()
    viewer.resize(900, 650)
    viewer.show()
    try:
        deadline = time.monotonic() + 30
        from PyQt6.QtTest import QTest
        while not viewer._ready and viewer.available and time.monotonic() < deadline:
            qt_app.processEvents()
            QTest.qWait(20)
        assert viewer._ready, viewer.status.text()
        big_receptor = "REMARK " + "x" * (2200000) + "\n" + PDB
        # A supplied contact exercises frontend lines/labels separately from the
        # chemistry detector's controlled geometry fixtures.
        contact = dict(kind="hbond", start=dict(x=1, y=2, z=3), end=dict(x=3.2, y=3.2, z=3),
                       distance=2.51, residue="ALA A:1", protein_serials=[1], detail="H-bond candidate")
        viewer.set_structure(ViewerStructure(big_receptor, ligand, "Target", "ethanol", "state", "run", "pose", 1, -5,
                                             [contact], "Geometric predictions"))
        wait_for(qt_app, lambda: viewer.export.isEnabled() or viewer.status.text().startswith("Viewer error"))
        assert viewer.export.isEnabled(), viewer.status.text()
        QTest.qWait(300)  # Let Qt and Chromium finish the label/viewport resize.
        width, height = viewer.web.width(), viewer.web.height()
        target = tmp_path / "native.png"
        viewer._export_path = str(target)
        viewer._send("export")
        wait_for(qt_app, target.exists)
        image = QImage(str(target))
        assert not image.isNull()
        assert (image.width(), image.height()) == (width * 2, height * 2)
        colors = {image.pixelColor(x, y).name() for x in range(0, image.width(), 8)
                  for y in range(0, image.height(), 8)}
        assert len(colors) > 10, "PNG should contain the molecule, not a blank canvas"
        viewer.interaction_toggles["hbond"].setChecked(False)
        without = tmp_path / "native-without-contacts.png"
        viewer._export_path = str(without)
        viewer._send("export")
        wait_for(qt_app, without.exists)
        assert QImage(str(without)) != image, "Contact lines/labels must change the exported molecular canvas"
    finally:
        viewer.close()
