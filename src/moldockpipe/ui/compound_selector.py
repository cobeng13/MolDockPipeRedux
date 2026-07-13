from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QDialog, QDialogButtonBox, QHBoxLayout, QLabel, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout

from ..project import ProjectRepository


class CompoundSelectorDialog(QDialog):
    def __init__(self, compounds: list[object], selected: set[str], repository: ProjectRepository, parent: QDialog | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Select Docked Compounds")
        self.resize(900, 650)
        self.table = QTableWidget(len(compounds), 4)
        self.table.setHorizontalHeaderLabels(["", "ID", "Vina score", "Exported Status"])
        self.table.setSortingEnabled(False)
        for row_index, compound in enumerate(compounds):
            parent_id = str(compound["parent_id"])
            check = QTableWidgetItem()
            check.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
            check.setCheckState(Qt.CheckState.Unchecked)
            self.table.setItem(row_index, 1, QTableWidgetItem(parent_id))
            self.table.setItem(row_index, 2, QTableWidgetItem(f"{float(compound['best_affinity']):.3f}"))
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
        for row in range(self.table.rowCount()):
            item = self.table.item(row, 0)
            if item and item.flags() & Qt.ItemFlag.ItemIsUserCheckable:
                item.setCheckState(state)

    def selected(self) -> list[str]:
        selected = []
        for row in range(self.table.rowCount()):
            item = self.table.item(row, 0)
            if item and item.checkState() == Qt.CheckState.Checked:
                selected.append(self.table.item(row, 1).text())
        return selected
