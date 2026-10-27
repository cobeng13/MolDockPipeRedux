import time

import pytest

from test_interaction_analysis import project, service
from test_results_viewer import qt_app


def wait(qt_app, predicate):
    deadline = time.monotonic() + 10
    while not predicate():
        qt_app.processEvents()
        time.sleep(.01)
        assert time.monotonic() < deadline


def test_report_dialog_scope_numeric_sort_and_download(project, qt_app, monkeypatch, tmp_path):
    from moldockpipe.ui.interaction_dialog import InteractionDialog
    repo, _ = project
    service(repo).run(profile_id="r")
    dialog = InteractionDialog(repo, "r", "run")
    wait(qt_app, lambda: dialog.reader is None)
    assert dialog.scope.count() == 2
    assert dialog.table.rowCount() == 2
    dialog.sort(2)  # affinity ascending, sorted by SQLite numeric values
    wait(qt_app, lambda: dialog.reader is None)
    assert dialog.table.item(0, 2).text() == "-7"
    destination = tmp_path / "report.csv"
    monkeypatch.setattr("moldockpipe.ui.interaction_dialog.QFileDialog.getSaveFileName", lambda *a: (str(destination), ""))
    dialog.download("interaction_summary.csv")
    assert "reference_retention_percent" in destination.read_text(encoding="utf-8-sig")
    dialog.done(0)


def test_worker_progress_cancel_and_close(project, qt_app, monkeypatch):
    from moldockpipe.ui.interaction_dialog import InteractionDialog
    from moldockpipe.pipeline import PipelineRunner
    repo, _ = project
    def run(runner, **kwargs):
        while not runner.cancelled():
            time.sleep(.01)
        raise InterruptedError("Interaction analysis cancelled")
    monkeypatch.setattr(PipelineRunner, "run_interactions", run)
    dialog = InteractionDialog(repo, "r", "run")
    wait(qt_app, lambda: dialog.reader is None)
    dialog.start()
    assert dialog.cancel.isEnabled()
    dialog.done(0)
    wait(qt_app, lambda: dialog.worker is None)
    assert dialog._closing


def test_opt_in_setting_defaults_off(project, qt_app):
    from moldockpipe.ui.settings_dialog import SettingsDialog
    repo, _ = project
    dialog = SettingsDialog(repo.get_settings(), repository=repo)
    assert not dialog.fields["postdock.analyze_interactions"].isChecked()
    dialog.fields["postdock.analyze_interactions"].setChecked(True)
    assert dialog.values()["postdock"]["analyze_interactions"] is True
    dialog.close()


def test_load_crystallographic_reference_action(project, qt_app, monkeypatch):
    from moldockpipe.ui.interaction_dialog import InteractionDialog, QMessageBox
    repo, _ = project
    dialog = InteractionDialog(repo, "r", "run")
    wait(qt_app, lambda: dialog.reader is None)
    path = repo.root / "inputs/receptors/r/reference.sdf"
    monkeypatch.setattr("moldockpipe.ui.interaction_dialog.QFileDialog.getOpenFileName", lambda *a, **kw: (str(path), ""))
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **kw: QMessageBox.StandardButton.Yes)
    dialog.attach_reference()
    wait(qt_app, lambda: dialog.reference_worker is None)
    assert repo.get_receptor_profiles()[0]["interaction_reference"]["source_kind"] == "manual_crystallographic_sdf"
    assert "source coordinates verified" in dialog.reference_label.text()
    assert "reference attached" in dialog.status.text()
    dialog.done(0)


def test_import_package_rebases_reference_paths_without_validation_history(project, qt_app, monkeypatch, tmp_path):
    import json
    import zipfile
    from moldockpipe.project import ProjectRepository
    from moldockpipe.ui.receptor_manager import ReceptorManagerDialog
    from moldockpipe.services.interaction_analysis import attach_crystallographic_reference, resolve_inputs
    source, _ = project
    attach_crystallographic_reference(source, "r", source.root / "inputs/receptors/r/reference.sdf")
    profile = source.get_receptor_profiles()[0]
    package = tmp_path / "receptor.zip"
    with zipfile.ZipFile(package, "w") as archive:
        archive.writestr("profile.json", json.dumps({k: v for k, v in profile.items() if k != "id"}))
        archive.write(source.root / profile["receptor"], "receptor.pdbqt")
        directory = source.root / "inputs/receptors/r"
        for path in directory.rglob("*"):
            if path.is_file() and path.name != "receptor.pdbqt":
                archive.write(path, "preparation/" + path.relative_to(directory).as_posix())
    target = ProjectRepository.create(tmp_path / "imported")
    monkeypatch.setattr("moldockpipe.ui.receptor_manager.QFileDialog.getOpenFileName", lambda *a, **kw: (str(package), ""))
    dialog = ReceptorManagerDialog(target)
    dialog._import_profile()
    dialog._save_locked()
    imported = target.get_receptor_profiles()[0]
    assert imported["id"] != "r"
    for ref in (imported["reference_ligand"], imported["interaction_reference"]):
        assert imported["id"] in ref["sdf"]
        assert (target.root / ref["sdf"]).is_file()
    # Source-evidence dictionary keys must be rebased too.
    assert all(imported["id"] in name for name in imported["interaction_reference"]["evidence_sha256"])
    assert resolve_inputs(target, imported, [])[1].GetNumAtoms() == 2
    imported.pop("interaction_reference")
    assert resolve_inputs(target, imported, [])[1].GetNumAtoms() == 2
    with target.connection() as conn:
        assert conn.execute("SELECT count(*) FROM redocking_runs").fetchone()[0] == 0
