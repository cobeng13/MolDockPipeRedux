from __future__ import annotations

import math
from PyQt6.QtCore import QAbstractTableModel, QModelIndex, Qt
from PyQt6.QtWidgets import (QAbstractItemView, QDialog, QDialogButtonBox, QDoubleSpinBox,
    QHBoxLayout, QLabel, QMessageBox, QPushButton, QTableView, QVBoxLayout)


class CompoundModel(QAbstractTableModel):
    headers = ("", "ID", "Best docking score", "Exported Status")

    def __init__(self, rows, selected=(), parent=None):
        super().__init__(parent)
        self.rows = [dict(row) for row in rows]
        self.checked = set(selected) & {str(row['parent_id']) for row in self.rows if not row.get('exported')}

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.rows)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else 4

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if role == Qt.ItemDataRole.DisplayRole and orientation == Qt.Orientation.Horizontal:
            return self.headers[section]

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        row = self.rows[index.row()]
        if index.column() == 0 and role == Qt.ItemDataRole.CheckStateRole:
            return Qt.CheckState.Checked if str(row['parent_id']) in self.checked else Qt.CheckState.Unchecked
        if role == Qt.ItemDataRole.DisplayRole:
            return ('', str(row['parent_id']), f"{float(row['best_affinity']):.3f}",
                    'Exported' if row.get('exported') else 'Not exported')[index.column()]
        if role == Qt.ItemDataRole.UserRole:
            return float(row['best_affinity'])

    def flags(self, index):
        flags = Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable
        if index.isValid() and index.column() == 0 and not self.rows[index.row()].get('exported'):
            flags |= Qt.ItemFlag.ItemIsUserCheckable
        return flags

    def setData(self, index, value, role=Qt.ItemDataRole.EditRole):
        if role != Qt.ItemDataRole.CheckStateRole or not self.flags(index) & Qt.ItemFlag.ItemIsUserCheckable:
            return False
        parent_id = str(self.rows[index.row()]['parent_id'])
        if value == Qt.CheckState.Checked or value == Qt.CheckState.Checked.value:
            self.checked.add(parent_id)
        else:
            self.checked.discard(parent_id)
        self.dataChanged.emit(index, index, [role])
        return True

    def sort(self, column, order=Qt.SortOrder.AscendingOrder):
        self.beginResetModel()
        def key(row):
            return (str(row['parent_id']) in self.checked, str(row['parent_id']),
                    float(row['best_affinity']), bool(row.get('exported')))[column]
        self.rows.sort(key=key, reverse=order == Qt.SortOrder.DescendingOrder)
        self.endResetModel()

    def select_where(self, predicate):
        self.checked = {str(row['parent_id']) for row in self.rows if not row.get('exported') and predicate(row)}
        if self.rows:
            self.dataChanged.emit(self.index(0, 0), self.index(len(self.rows)-1, 0), [Qt.ItemDataRole.CheckStateRole])


class CompoundSelectorDialog(QDialog):
    def __init__(self, compounds, selected, repository, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select Docked Compounds")
        self.resize(900, 650)
        self.model = CompoundModel(compounds, selected, self)
        self.table = QTableView()
        self.table.setModel(self.model)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSortingEnabled(True)
        self.table.setColumnWidth(1, 240)
        self.table.setColumnWidth(2, 180)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Select by best docking score across enabled receptors (kcal/mol, inclusive)."))
        controls = QHBoxLayout()
        self.score_start, self.score_end = QDoubleSpinBox(), QDoubleSpinBox()
        for field in (self.score_start, self.score_end):
            field.setDecimals(3); field.setRange(-1000000, 1000000); field.setSingleStep(0.5)
        scores = [float(row['best_affinity']) for row in self.model.rows]
        if scores:
            self.score_start.setValue(math.floor(min(scores)*1000)/1000)
            self.score_end.setValue(math.ceil(max(scores)*1000)/1000)
        controls.addWidget(QLabel("Start")); controls.addWidget(self.score_start)
        controls.addWidget(QLabel("End")); controls.addWidget(self.score_end)
        self.select_range_button = QPushButton("Select score range")
        self.select_range_button.clicked.connect(self._select_score_range)
        self.select_range_button.setEnabled(bool(compounds))
        controls.addWidget(self.select_range_button)
        layout.addLayout(controls)
        self.selection_count = QLabel()
        layout.addWidget(self.selection_count)
        layout.addWidget(self.table)
        self.model.dataChanged.connect(self._update_selection_count)
        self._update_selection_count()
        marks = QHBoxLayout()
        for label, checked in (("Mark All", True), ("Unmark All", False)):
            button = QPushButton(label)
            button.clicked.connect(lambda _, value=checked: self._set_all_checked(value))
            marks.addWidget(button)
        layout.addLayout(marks)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept); buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _update_selection_count(self, *_):
        self.selection_count.setText(f"{len(self.model.checked)} compounds selected for export")

    def _set_all_checked(self, checked):
        self.model.select_where(lambda row: checked)
        self._update_selection_count()

    def _select_score_range(self):
        start, end = self.score_start.value(), self.score_end.value()
        if start > end:
            QMessageBox.warning(self, "Invalid score range", "Start must be less than or equal to End (for example, -10 to -7).")
            return
        self.model.select_where(lambda row: start <= float(row['best_affinity']) <= end)
        self._update_selection_count()

    def selected(self):
        return sorted(self.model.checked)
