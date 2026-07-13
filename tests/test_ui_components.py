from __future__ import annotations

import pytest

try:
    from moldockpipe.ui.compound_selector import CompoundSelectorDialog
    from moldockpipe.ui.progress import CheckpointProgress
    from moldockpipe.ui.settings_dialog import SettingsDialog
    from moldockpipe.ui.workers import PipelineWorker
    import moldockpipe.ui.workers as workers
except ImportError as exc:
    pytest.skip(f"Qt UI modules unavailable: {exc}", allow_module_level=True)


def test_ui_components_import_and_expose_expected_contracts() -> None:
    assert SettingsDialog.workflow_reset is not None
    assert CompoundSelectorDialog.selected
    assert CheckpointProgress.STAGES == ("screening", "molscrub", "meeko", "vina", "postdock")


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
    worker = PipelineWorker(object(), ["meeko"])
    worker.finished.connect(summaries.append)
    worker.run()

    assert calls == ["meeko"]
    assert summaries == [{"meeko": (1, 0)}]
