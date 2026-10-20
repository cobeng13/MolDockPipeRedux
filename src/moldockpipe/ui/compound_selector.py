from __future__ import annotations

import math

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QAbstractItemView, QDialog, QDialogButtonBox, QDoubleSpinBox, QHBoxLayout, QLabel, QMessageBox, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout

from ..project import ProjectRepository


class CompoundSelectorDialog(QDialog):
    def __init__(self, compounds: list[object], selected: set[str], repository: ProjectRepository, parent: QDialog | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Select Docked Compounds")
        self.resize(900, 650)
        self.table = QTableWidget(len(compounds), 4)
        self.table.setHorizontalHeaderLabels(["", "ID", "Best docking score", "Exported Status"])
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSortingEnabled(False)
        for row_index, compound in enumerate(compounds):
            parent_id = str(compound["parent_id"])
            check = QTableWidgetItem()
            check.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
            check.setCheckState(Qt.CheckState.Unchecked)
            self.table.setItem(row_index, 1, QTableWidgetItem(parent_id))
            score = float(compound['best_affinity'])
            score_item = QTableWidgetItem(f"{score:.3f}")
            score_item.setData(Qt.ItemDataRole.UserRole, score)
            self.table.setItem(row_index, 2, score_item)
            exported = any((repository.root / folder / parent_id).exists() for folder in ("For_PostDocking/SDF", "For_PostDocking/PDBQTs"))
            status = QTableWidgetItem("Exported" if exported else "Not exported")
            if exported:
                check.setFlags(Qt.ItemFlag.ItemIsEnabled)
                check.setCheckState(Qt.CheckState.Unchecked)
            self.table.setItem(row_index, 0, check)
            self.table.setItem(row_index, 3, status)
        self.table.resizeColumnsToContents()
        self.table.setSortingEnabled(True)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Select by best docking score across enabled receptors (kcal/mol, inclusive)."))
        range_controls = QHBoxLayout()
        self.score_start = QDoubleSpinBox()
        self.score_end = QDoubleSpinBox()
        for field in (self.score_start, self.score_end):
            field.setDecimals(3)
            field.setRange(-1000000, 1000000)
            field.setSingleStep(0.5)
        scores = [float(compound['best_affinity']) for compound in compounds]
        if scores:
            # Round outward so the initial range includes scores with more precision.
            self.score_start.setValue(math.floor(min(scores) * 1000) / 1000)
            self.score_end.setValue(math.ceil(max(scores) * 1000) / 1000)
        range_controls.addWidget(QLabel("Start"))
        range_controls.addWidget(self.score_start)
        range_controls.addWidget(QLabel("End"))
        range_controls.addWidget(self.score_end)
        self.select_range_button = QPushButton("Select score range")
        self.select_range_button.clicked.connect(self._select_score_range)
        self.select_range_button.setEnabled(bool(compounds))
        range_controls.addWidget(self.select_range_button)
        layout.addLayout(range_controls)
        self.selection_count = QLabel()
        layout.addWidget(self.selection_count)
        self.table.itemChanged.connect(self._update_selection_count)
        self._update_selection_count()
        layout.addWidget(self.table)
        if not compounds:
            layout.addWidget(QLabel("No successfully docked compounds are available."))
        selection_buttons = QHBoxLayout()
        mark_all = QPushButton("Mark All")
        mark_all.clicked.connect(lambda: self._set_all_checked(True))
        unmark_all = QPushButton("Unmark All")
        unmark_all.clicked.connect(lambda: self._set_all_checked(False))
        selection_buttons.addWidget(mark_all)
        selection_buttons.addWidget(unmark_all)
        selection_buttons.addStretch()
        layout.addLayout(selection_buttons)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _set_all_checked(self, checked: bool) -> None:
        state = Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
        sorting = self.table.isSortingEnabled()
        self.table.setSortingEnabled(False)
        self.table.blockSignals(True)
        try:
            for row in range(self.table.rowCount()):
                item = self.table.item(row, 0)
                if item and item.flags() & Qt.ItemFlag.ItemIsUserCheckable:
                    item.setCheckState(state)
        finally:
            self.table.setSortingEnabled(sorting)
            self.table.blockSignals(False)
        self._update_selection_count()

    def _select_score_range(self) -> None:
        start, end = self.score_start.value(), self.score_end.value()
        if start > end:
            QMessageBox.warning(self, "Invalid score range", "Start must be less than or equal to End (for example, -10 to -7).")
            return
        sorting = self.table.isSortingEnabled()
        self.table.setSortingEnabled(False)
        self.table.blockSignals(True)
        try:
            for row in range(self.table.rowCount()):
                check = self.table.item(row, 0)
                if check.flags() & Qt.ItemFlag.ItemIsUserCheckable:
                    score = self.table.item(row, 2).data(Qt.ItemDataRole.UserRole)
                    check.setCheckState(Qt.CheckState.Checked if start <= score <= end else Qt.CheckState.Unchecked)
        finally:
            self.table.setSortingEnabled(sorting)
            self.table.blockSignals(False)
        self._update_selection_count()

    def _update_selection_count(self) -> None:
        self.selection_count.setText(f"{len(self.selected())} compounds selected for export")

    def selected(self) -> list[str]:
        selected = []
        for row in range(self.table.rowCount()):
            item = self.table.item(row, 0)
            if item and item.checkState() == Qt.CheckState.Checked:
                selected.append(self.table.item(row, 1).text())
        return selected
