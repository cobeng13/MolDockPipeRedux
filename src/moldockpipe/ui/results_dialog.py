from __future__ import annotations

import json

from PyQt6.QtCore import QThread, Qt, pyqtSignal
from PyQt6.QtWidgets import QAbstractItemView, QDialog, QDialogButtonBox, QHBoxLayout, QLabel, QPushButton, QSplitter, QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget

from ..project import ProjectRepository
from ..services.structure_viewer import StructureViewerService
from .molecular_viewer import MolecularViewer
from .workers import TaskWorker


def load_profile_page(repository, profile_id, offset, limit=200, cancelled=lambda: False):
    with repository.connection() as conn:
        conn.set_progress_handler(lambda: int(cancelled()), 5000)
        total = conn.execute("SELECT COUNT(*) FROM parent_ligands").fetchone()[0]
        rows = conn.execute("""WITH ranked AS (
            SELECT s.parent_id,s.state_id,p.affinity,p.mode_index,d.command_json,d.run_id,p.pose_id,
                ROW_NUMBER() OVER (PARTITION BY s.parent_id ORDER BY p.affinity,s.state_id,p.mode_index,d.run_id,p.pose_id) rank
            FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id
            JOIN docking_poses p ON p.run_id=d.run_id
            WHERE d.receptor_profile_id=? AND d.is_current=1 AND d.status='completed' AND s.active=1),
        latest AS (
            SELECT s.parent_id,d.status,d.command_json,COALESCE(d.reason,'') reason,
                ROW_NUMBER() OVER (PARTITION BY s.parent_id ORDER BY d.started_at DESC,d.run_id DESC) rank
            FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id
            WHERE d.receptor_profile_id=? AND d.is_current=1 AND s.active=1)
        SELECT l.parent_id,b.state_id,b.affinity,b.mode_index,b.run_id,b.pose_id,
            CASE WHEN b.run_id IS NOT NULL THEN 'completed' ELSE COALESCE(r.status,'pending') END status,
            CASE WHEN b.run_id IS NOT NULL THEN '' ELSE COALESCE(r.reason,'') END reason,
            COALESCE(b.command_json,r.command_json,'{}') command_json
        FROM parent_ligands l LEFT JOIN ranked b ON b.parent_id=l.parent_id AND b.rank=1
        LEFT JOIN latest r ON r.parent_id=l.parent_id AND r.rank=1
        ORDER BY b.affinity IS NULL,b.affinity,l.parent_id LIMIT ? OFFSET ?""",
            (profile_id, profile_id, limit, offset)).fetchall()
        conn.set_progress_handler(None, 0)
    return total, [dict(row) for row in rows]


class StructureLoadWorker(QThread):
    loaded = pyqtSignal(int, object)
    failed = pyqtSignal(int, str)

    def __init__(self, service, generation, run_id, pose_id, parent=None):
        super().__init__(parent)
        self.service, self.generation = service, generation
        self.run_id, self.pose_id = run_id, pose_id

    def run(self):
        try:
            result = self.service.load(self.run_id, self.pose_id, cancelled=self.isInterruptionRequested)
            if not self.isInterruptionRequested():
                self.loaded.emit(self.generation, result)
        except Exception as exc:
            if not self.isInterruptionRequested():
                self.failed.emit(self.generation, str(exc))


class DockingResultsDialog(QDialog):
    def __init__(self, repository: ProjectRepository, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.repository = repository
        self.service = StructureViewerService(repository)
        self._generation = 0
        self._worker = None
        self._table_worker = None
        self._table_pending = None
        self._pending = None
        self._closing = False
        self._run_id = None
        self.setWindowTitle("Multidock Results")
        self.resize(1350, 700)
        self.tabs = tabs = QTabWidget()
        self.viewer = MolecularViewer(self)
        profiles = repository.get_receptor_profiles(include_archived=True)
        for profile in profiles:
            tabs.addTab(self._profile_table(str(profile["id"])), str(profile.get("name", profile["id"])))
        if not profiles:
            tabs.addTab(QLabel("No receptor profiles configured."), "Results")
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(tabs)
        splitter.addWidget(self.viewer)
        splitter.setSizes([650, 700])
        splitter.setChildrenCollapsible(False)
        tabs.currentChanged.connect(self._tab_changed)
        self.viewer.pose.currentIndexChanged.connect(self._pose_changed)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close); buttons.rejected.connect(self.reject)
        self.previous_page = QPushButton("Previous 200")
        self.next_page = QPushButton("Next 200")
        self.page_label = QLabel()
        self.previous_page.clicked.connect(lambda: self._change_page(-200))
        self.next_page.clicked.connect(lambda: self._change_page(200))
        pages = QHBoxLayout()
        pages.addWidget(self.previous_page); pages.addWidget(self.page_label); pages.addWidget(self.next_page)
        self.interaction_action = QPushButton("Analyze Interactions")
        self.interaction_action.clicked.connect(self._open_interactions)
        pages.addWidget(self.interaction_action)
        layout = QVBoxLayout(self); layout.addWidget(splitter); layout.addLayout(pages); layout.addWidget(buttons)
        self._tab_changed()

    def _open_interactions(self):
        table = self.tabs.currentWidget()
        if not isinstance(table, QTableWidget):
            return
        from .interaction_dialog import InteractionDialog
        InteractionDialog(self.repository, table._profile_id, self._run_id, self).exec()

    def _tab_changed(self, *_):
        table = self.tabs.currentWidget()
        if isinstance(table, QTableWidget):
            if not table._loaded:
                self._load_page(table)
                return
            self._update_page_controls(table)
            if table.currentRow() < 0 and table.rowCount():
                table.selectRow(0)
            else:
                self._selection_changed(table)
        else:
            self._invalidate("No receptor profiles configured.")

    def _invalidate(self, message):
        self._generation += 1
        self._pending = None
        if self._worker:
            self._worker.requestInterruption()
        self.viewer.clear(message)

    def _selection_changed(self, table):
        if self.tabs.currentWidget() is not table:
            return
        self._invalidate("Select a completed docking result.")
        self._run_id = None
        self.viewer.pose.blockSignals(True)
        self.viewer.pose.clear()
        self.viewer.pose.setEnabled(False)
        item = table.item(table.currentRow(), 0)
        record = item.data(Qt.ItemDataRole.UserRole) if item else None
        try:
            if not record:
                if item:
                    status = table.item(table.currentRow(), 4).text()
                    reason = table.item(table.currentRow(), 5).text()
                    self.viewer.clear(f"No docking pose available ({status}). {reason}")
                return
            self._run_id = record["run_id"]
            for pose in self.service.poses(self._run_id):
                self.viewer.pose.addItem(f"Mode {pose['mode_index']} · {pose['affinity']:g} kcal/mol", pose["pose_id"])
            index = self.viewer.pose.findData(record["pose_id"])
            if index < 0:
                raise ValueError("The selected docking pose is no longer available.")
            self.viewer.pose.setCurrentIndex(index)
            self.viewer.pose.setEnabled(self.viewer.available)
        except Exception as exc:
            self._run_id = None
            self.viewer.clear(f"Unable to load result: {exc}")
        finally:
            self.viewer.pose.blockSignals(False)
        if self._run_id:
            self._pose_changed()

    def _pose_changed(self, *_):
        self._invalidate("Loading receptor and ligand…")
        pose_id = self.viewer.pose.currentData()
        if self._run_id and pose_id and self.viewer.available:
            self._pending = (self._generation, self._run_id, pose_id)
            self._start_pending()

    def _start_pending(self):
        if self._worker is not None or self._pending is None or self._closing:
            return
        generation, run_id, pose_id = self._pending
        self._pending = None
        self._worker = StructureLoadWorker(self.service, generation, run_id, pose_id, self)
        self._worker.loaded.connect(self._loaded)
        self._worker.failed.connect(self._failed)
        self._worker.finished.connect(self._worker_finished)
        self._worker.start()

    def _loaded(self, generation, structure):
        if generation == self._generation and not self._closing and self.viewer.available:
            self.viewer.set_structure(structure)

    def _failed(self, generation, message):
        if generation == self._generation and not self._closing:
            self.viewer.clear(f"Unable to load pose: {message}")

    def _worker_finished(self):
        worker, self._worker = self._worker, None
        worker.deleteLater()
        if self._closing:
            self._finish_close()
        else:
            self._start_pending()

    def done(self, result):
        self._closing = True
        self._close_result = result
        self._invalidate("Closing viewer…")
        self._table_pending = None
        if self._table_worker is not None:
            self._table_worker.requestInterruption()
        if self._worker is None and self._table_worker is None:
            super().done(result)
        else:
            self.setEnabled(False)

    def closeEvent(self, event):
        self.done(QDialog.DialogCode.Rejected)
        if self._worker is not None or self._table_worker is not None:
            event.ignore()
        else:
            event.accept()

    def _finish_close(self):
        if self._worker is None and self._table_worker is None:
            super().done(self._close_result)

    def _profile_table(self, profile_id: str) -> QTableWidget:
        headers = ["Ligand", "Best state", "Affinity", "Mode", "Status", "Reason", "Protocol"]
        table = QTableWidget(0, len(headers))
        table.setHorizontalHeaderLabels(headers)
        table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        table._profile_id, table._offset, table._total, table._loaded = profile_id, 0, 0, False
        table.itemSelectionChanged.connect(lambda: self._selection_changed(table))
        return table

    def _update_page_controls(self, table):
        self.previous_page.setEnabled(table._offset > 0)
        self.next_page.setEnabled(table._offset + 200 < table._total)
        first = table._offset + 1 if table._total else 0
        self.page_label.setText(f"{first}-{min(table._offset + 200, table._total)} of {table._total}")

    def _change_page(self, amount):
        table = self.tabs.currentWidget()
        if isinstance(table, QTableWidget):
            table._offset = max(0, table._offset + amount)
            table._loaded = False
            self._load_page(table)

    def _load_page(self, table):
        self._invalidate("Loading results...")
        self._run_id = None
        self.previous_page.setEnabled(False); self.next_page.setEnabled(False)
        self.page_label.setText("Loading...")
        self._table_pending = (table, table._offset)
        if self._table_worker is not None:
            self._table_worker.requestInterruption()
        self._start_table_pending()

    def _start_table_pending(self):
        if self._table_worker is not None or self._table_pending is None or self._closing:
            return
        table, offset = self._table_pending
        self._table_pending = None
        profile_id = table._profile_id
        worker = TaskWorker(lambda: load_profile_page(self.repository, profile_id, offset,
                                                      cancelled=worker.isInterruptionRequested), self)
        self._table_worker = worker
        worker.succeeded.connect(lambda result: self._page_ready(table, offset, result))
        worker.failed.connect(lambda message: self._page_failed(table, offset, message))
        worker.finished.connect(self._table_finished)
        worker.start()

    def _page_failed(self, table, offset, message):
        if not self._closing and self.tabs.currentWidget() is table and table._offset == offset and self._table_pending is None:
            self.page_label.setText("Unable to load results")
            self.previous_page.setEnabled(offset > 0)
            self.viewer.clear(f"Unable to load results: {message}")

    def _page_ready(self, table, offset, result):
        if self._closing or self.tabs.currentWidget() is not table or table._offset != offset or self._table_pending is not None:
            return
        total, rows = result
        table.blockSignals(True)
        table.clearSelection()
        table.setCurrentCell(-1, -1)
        table.setRowCount(len(rows))
        for i, row in enumerate(rows):
            values = [row['parent_id'], row['state_id'] or '', row['affinity'] if row['affinity'] is not None else '',
                      row['mode_index'] or '', row['status'], row['reason'],
                      json.loads(row['command_json']).get('settings', {}).get('protocol', 'vina')]
            for column, value in enumerate(values):
                table.setItem(i, column, QTableWidgetItem(str(value)))
            if row['pose_id']:
                table.item(i, 0).setData(Qt.ItemDataRole.UserRole,
                    {key: row[key] for key in ('run_id', 'pose_id', 'state_id')})
        table._total, table._loaded = total, True
        table.blockSignals(False)
        self._update_page_controls(table)
        if rows:
            table.selectRow(0)
        else:
            self.viewer.clear("No docking results available.")

    def _table_finished(self):
        worker, self._table_worker = self._table_worker, None
        worker.deleteLater()
        if self._closing:
            self._finish_close()
        else:
            self._start_table_pending()
