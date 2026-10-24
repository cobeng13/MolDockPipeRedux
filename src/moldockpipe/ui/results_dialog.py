from __future__ import annotations

import json

from PyQt6.QtCore import QThread, Qt, pyqtSignal
from PyQt6.QtWidgets import QAbstractItemView, QDialog, QDialogButtonBox, QLabel, QSplitter, QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget

from ..project import ProjectRepository
from ..services.structure_viewer import StructureViewerService
from .molecular_viewer import MolecularViewer


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
        layout = QVBoxLayout(self); layout.addWidget(splitter); layout.addWidget(buttons)
        self._tab_changed()

    def _tab_changed(self, *_):
        table = self.tabs.currentWidget()
        if isinstance(table, QTableWidget):
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
            super().done(self._close_result)
        else:
            self._start_pending()

    def done(self, result):
        self._closing = True
        self._close_result = result
        self._invalidate("Closing viewer…")
        if self._worker is None:
            super().done(result)
        else:
            self.setEnabled(False)

    def closeEvent(self, event):
        self.done(QDialog.DialogCode.Rejected)
        if self._worker is not None:
            event.ignore()
        else:
            event.accept()

    def _profile_table(self, profile_id: str) -> QTableWidget:
        headers = ["Ligand", "Best state", "Affinity", "Mode", "Status", "Reason", "Protocol"]
        table = QTableWidget(0, len(headers)); table.setHorizontalHeaderLabels(headers)
        table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        with self.repository.connection() as conn:
            parents = conn.execute("SELECT parent_id FROM parent_ligands ORDER BY parent_id").fetchall()
            best = {row["parent_id"]: row for row in conn.execute("""WITH ranked AS (
                    SELECT s.parent_id, s.state_id, p.affinity, p.mode_index, d.command_json, d.run_id, p.pose_id,
                           ROW_NUMBER() OVER (PARTITION BY s.parent_id ORDER BY p.affinity, s.state_id, p.mode_index, d.run_id, p.pose_id) rank
                    FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id
                    JOIN docking_poses p ON p.run_id=d.run_id
                    WHERE d.receptor_profile_id=? AND d.is_current=1 AND d.status='completed' AND s.active=1)
                SELECT * FROM ranked WHERE rank=1""", (profile_id,)).fetchall()}
            latest = {row["parent_id"]: row for row in conn.execute("""WITH ranked AS (
                    SELECT s.parent_id, d.status, d.command_json, COALESCE(d.reason, '') reason,
                           ROW_NUMBER() OVER (PARTITION BY s.parent_id ORDER BY d.started_at DESC, d.run_id DESC) rank
                    FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id
                    WHERE d.receptor_profile_id=? AND d.is_current=1 AND s.active=1)
                SELECT * FROM ranked WHERE rank=1""", (profile_id,)).fetchall()}
        parents = sorted(parents, key=lambda parent: (
            parent["parent_id"] not in best,
            best[parent["parent_id"]]["affinity"] if parent["parent_id"] in best else 0,
            parent["parent_id"],
        ))
        for parent in parents:
            parent_id = parent["parent_id"]
            score = best.get(parent_id); run = latest.get(parent_id)
            row = table.rowCount(); table.insertRow(row)
            values = [parent_id, score["state_id"] if score else "", score["affinity"] if score else "",
                      score["mode_index"] if score else "", "completed" if score else run["status"] if run else "pending",
                      "" if score else run["reason"] if run else ""]
            source = score or run
            values.append(json.loads(source["command_json"]).get("settings", {}).get("protocol", "vina") if source else "")
            for column, value in enumerate(values):
                table.setItem(row, column, QTableWidgetItem(str(value)))
            if score:
                table.item(row, 0).setData(Qt.ItemDataRole.UserRole, {
                    "run_id": score["run_id"], "pose_id": score["pose_id"], "state_id": score["state_id"],
                })
        table.itemSelectionChanged.connect(lambda: self._selection_changed(table))
        return table
