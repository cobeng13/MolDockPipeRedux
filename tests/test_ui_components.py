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


def test_pipeline_worker_routes_requested_stage(monkeypatch) -> None:
    calls: list[str] = []

    class FakeRunner:
        def __init__(self, repository, progress=None) -> None:
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
    worker = PipelineWorker(SimpleNamespace(import_project_input=lambda: calls.append("import")), ["meeko"])
    worker.finished.connect(summaries.append)
    worker.run()

    assert calls == ["import", "meeko"]
    assert summaries == [{"meeko": (1, 0)}]


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
