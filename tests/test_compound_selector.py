from __future__ import annotations

import pytest

pytest.importorskip("PyQt6")
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QApplication, QMessageBox

from moldockpipe.project import ProjectRepository
from moldockpipe.ui.compound_selector import CompoundSelectorDialog


@pytest.fixture
def selector(tmp_path):
    app = QApplication.instance() or QApplication([])
    repo = ProjectRepository.create(tmp_path / "project")
    compounds = [
        {"parent_id": name, "best_affinity": score}
        for name, score in [("lower", -10), ("middle", -8), ("upper", -7),
                            ("outside", -6.9996), ("exported", -9)]
    ]
    (repo.root / "For_PostDocking" / "SDF" / "exported").mkdir(parents=True)
    dialog = CompoundSelectorDialog(compounds, set(), repo)
    yield dialog
    dialog.close()
    app.processEvents()


def test_score_range_includes_endpoints_after_sort_and_replaces_selection(selector):
    selector._set_all_checked(True)
    selector.table.sortItems(1, Qt.SortOrder.DescendingOrder)
    selector.score_start.setValue(-10)
    selector.score_end.setValue(-7)
    selector.select_range_button.click()
    assert set(selector.selected()) == {"lower", "middle", "upper"}
    assert selector.selection_count.text() == "3 compounds selected for export"
    # The displayed rounded score of outside is -7.000, but its real value is outside.
    assert "outside" not in selector.selected()
    # Selecting again must remove the previously selected compounds.
    selector.score_start.setValue(-8)
    selector.score_end.setValue(-8)
    selector.select_range_button.click()
    assert selector.selected() == ["middle"]


def test_reversed_range_keeps_selection_and_empty_range_clears_it(selector, monkeypatch):
    selector._set_all_checked(True)
    before = set(selector.selected())
    warnings = []
    monkeypatch.setattr(QMessageBox, "warning", lambda *args: warnings.append(args[2]))
    selector.score_start.setValue(-7)
    selector.score_end.setValue(-10)
    selector.select_range_button.click()
    assert warnings and set(selector.selected()) == before
    selector.score_start.setValue(-20)
    selector.score_end.setValue(-15)
    selector.select_range_button.click()
    assert selector.selected() == []
    assert selector.selection_count.text() == "0 compounds selected for export"


def test_empty_selector_disables_range_selection(tmp_path):
    app = QApplication.instance() or QApplication([])
    repo = ProjectRepository.create(tmp_path / "empty")
    dialog = CompoundSelectorDialog([], set(), repo)
    assert not dialog.select_range_button.isEnabled()
    assert dialog.selected() == []
    dialog.close()
    app.processEvents()
