import os
from types import SimpleNamespace

import pytest

from moldockpipe.execution_lock import project_execution_lock
from moldockpipe.project import ProjectRepository

pytest.importorskip("PyQt6")
from moldockpipe.ui import settings_dialog as settings


@pytest.fixture
def project(tmp_path):
    repo = ProjectRepository.create(tmp_path / "project")
    source = repo.root / "inputs" / "input.csv"
    source.write_text("id,smiles\nlig1,CCO\n", encoding="utf-8")
    repo.import_ligands_csv(source)
    for name in ("artifacts", "logs", "exports", "For_PostDocking"):
        folder = repo.root / name
        folder.mkdir(exist_ok=True)
        (folder / "generated.txt").write_text("generated", encoding="utf-8")
    return repo


def parent_count(repo):
    with repo.connection() as conn:
        return conn.execute("SELECT COUNT(*) FROM parent_ligands").fetchone()[0]


def test_purge_removes_outputs_and_rows_but_preserves_inputs_and_settings(project):
    config = (project.root / "project.yml").read_bytes()
    settings.SettingsDialog._clear_workflow_files(SimpleNamespace(repository=project))
    assert parent_count(project) == 0
    for name in ("artifacts", "logs", "exports", "For_PostDocking"):
        assert list((project.root / name).iterdir()) == []
    assert (project.root / "inputs" / "input.csv").exists()
    assert (project.root / "project.yml").read_bytes() == config


def test_purge_retries_transient_windows_sharing_violation(project, monkeypatch):
    remove = settings.shutil.rmtree
    attempts = []

    def transient(path):
        attempts.append(path)
        if len(attempts) == 1:
            exc = PermissionError("file in use")
            exc.winerror = 32
            raise exc
        remove(path)

    monkeypatch.setattr(settings.shutil, "rmtree", transient)
    monkeypatch.setattr(settings.time, "sleep", lambda _: None)
    settings.SettingsDialog._clear_workflow_files(SimpleNamespace(repository=project))
    assert attempts[:2] == [project.root / "artifacts"] * 2
    assert parent_count(project) == 0


@pytest.mark.parametrize("winerror, expected_attempts", [(32, 3), (33, 3), (5, 1)])
def test_failed_cleanup_keeps_database_rows_and_releases_project_lock(
    project, monkeypatch, winerror, expected_attempts
):
    attempts = []

    def fail(path):
        attempts.append(path)
        exc = PermissionError("cannot remove file")
        exc.winerror = winerror
        raise exc

    monkeypatch.setattr(settings.shutil, "rmtree", fail)
    monkeypatch.setattr(settings.time, "sleep", lambda _: None)
    with pytest.raises(PermissionError):
        settings.SettingsDialog._clear_workflow_files(SimpleNamespace(repository=project))
    assert len(attempts) == expected_attempts
    assert parent_count(project) == 1
    with project_execution_lock(project.root):
        pass


def test_purge_reports_locked_file_without_success_or_reset(monkeypatch):
    warnings, completed, resets = [], [], []

    def warn(parent, title, message, *args):
        warnings.append((title, message))
        return settings.QMessageBox.StandardButton.Yes

    def fail():
        raise PermissionError(13, "file in use", "ligand.pdbqt")

    monkeypatch.setattr(settings.QMessageBox, "warning", warn)
    monkeypatch.setattr(settings.QMessageBox, "information", lambda *args: completed.append(args))
    fake = SimpleNamespace(repository=object(), _clear_workflow_files=fail,
                           workflow_reset=SimpleNamespace(emit=lambda: resets.append(True)))
    settings.SettingsDialog._purge_workflow(fake)
    assert warnings[-1][0] == "Purge incomplete"
    assert "ligand.pdbqt" in warnings[-1][1]
    assert completed == resets == []


@pytest.mark.skipif(os.name != "nt", reason="Windows denies deletion of an open file")
def test_purge_with_real_windows_file_lock_can_be_retried(project):
    ligand = project.root / "artifacts" / "ligand.pdbqt"
    ligand.write_text("ligand", encoding="utf-8")
    fake = SimpleNamespace(repository=project)
    with ligand.open("rb"):
        with pytest.raises(PermissionError) as error:
            settings.SettingsDialog._clear_workflow_files(fake)
        assert error.value.winerror == 32
        assert ligand.exists()
        assert parent_count(project) == 1
    settings.SettingsDialog._clear_workflow_files(fake)
    assert parent_count(project) == 0
    assert not ligand.exists()
