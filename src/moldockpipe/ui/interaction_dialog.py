"""Interaction reports in the existing Results workflow (no WebEngine required)."""
from __future__ import annotations

import json
import shutil

from PyQt6.QtCore import QThread, Qt, pyqtSignal
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QDialog, QFileDialog, QHBoxLayout,
    QLabel, QMessageBox, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout)

from ..pipeline import PipelineRunner
from .workers import TaskWorker


class InteractionWorker(QThread):
    progress = pyqtSignal(object)
    succeeded = pyqtSignal(object)
    failed = pyqtSignal(str)

    def __init__(self, repository, profile_id, run_id, force, parent=None):
        super().__init__(parent)
        self.repository, self.profile_id, self.run_id, self.force = repository, profile_id, run_id, force

    def run(self):
        try:
            result = PipelineRunner(self.repository, progress=self.progress.emit,
                cancelled=self.isInterruptionRequested).run_interactions(
                    profile_id=self.profile_id, run_id=self.run_id, force=self.force)
            self.succeeded.emit(result)
        except Exception as exc:
            self.failed.emit(str(exc))


class InteractionDialog(QDialog):
    COLUMNS = [("ligand_id", "Ligand"), ("pose_id", "Pose"), ("affinity", "Affinity"),
        ("reference_retention_percent", "Reference retained %"), ("jaccard_similarity", "Jaccard"),
        ("conserved_interaction_pairs", "Conserved pairs"), ("lost_reference_pairs", "Lost pairs"),
        ("new_interaction_pairs", "New pairs"), ("residue_retention_percent", "Residues retained %"),
        ("analysis_status", "Status"), ("error_message", "Error")]

    def __init__(self, repository, profile_id, selected_run=None, parent=None):
        super().__init__(parent)
        self.repository, self.profile_id, self.selected_run = repository, profile_id, selected_run
        self.worker = self.reader = None
        self.reference_worker = None
        self._closing = False
        self._read_pending = False
        self._last_error = None
        self._offset, self._sort, self._descending = 0, "reference_retention_percent", True
        self._output = None
        self.setWindowTitle("Protein–ligand interaction analysis")
        self.resize(1200, 600)
        layout = QVBoxLayout(self)
        self.scope = QComboBox()
        self.scope.addItem("All exported poses for this receptor", None)
        if selected_run:
            self.scope.addItem("Selected docking run", selected_run)
        self.force = QCheckBox("Create a new snapshot")
        self.analyze = QPushButton("Analyze Interactions")
        self.cancel = QPushButton("Cancel analysis")
        self.cancel.setEnabled(False)
        controls = QHBoxLayout()
        for widget in (self.scope, self.force, self.analyze, self.cancel):
            controls.addWidget(widget)
        layout.addLayout(controls)
        reference_controls = QHBoxLayout()
        self.reference_label = QLabel()
        self.reference_label.setWordWrap(True)
        self.load_reference = QPushButton("Load co-crystallized SDF…")
        self.load_reference.clicked.connect(self.attach_reference)
        reference_controls.addWidget(self.reference_label, 1)
        reference_controls.addWidget(self.load_reference)
        layout.addLayout(reference_controls)
        self.refresh_reference_label()
        layout.addWidget(QLabel("ProLIF interaction-pattern similarity; these metrics do not predict affinity or activity.\n"
            "The 3D viewer's geometric annotations are calculated separately."))
        self.status = QLabel("Loading saved analysis…")
        layout.addWidget(self.status)
        self.table = QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels([label for _, label in self.COLUMNS])
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSortIndicatorShown(True)
        self.table.horizontalHeader().sectionClicked.connect(self.sort)
        layout.addWidget(self.table)
        controls = QHBoxLayout()
        self.previous, self.next = QPushButton("Previous 200"), QPushButton("Next 200")
        self.previous.clicked.connect(lambda: self.page(-200))
        self.next.clicked.connect(lambda: self.page(200))
        controls.addWidget(self.previous); controls.addWidget(self.next)
        self.downloads = []
        for name in ("interactions_detailed.csv", "interaction_summary.csv", "interaction_comparison.csv"):
            button = QPushButton(f"Save {name}")
            button.clicked.connect(lambda checked=False, filename=name: self.download(filename))
            controls.addWidget(button)
            self.downloads.append(button)
        layout.addLayout(controls)
        self.analyze.clicked.connect(self.start)
        self.cancel.clicked.connect(lambda: self.worker.requestInterruption() if self.worker else None)
        self.scope.currentIndexChanged.connect(lambda: self.page(-self._offset))
        self.load()

    def refresh_reference_label(self):
        profile = next((p for p in self.repository.get_receptor_profiles(include_archived=True)
                        if str(p["id"]) == self.profile_id), {})
        record = profile.get("interaction_reference") or profile.get("reference_ligand") or {}
        origin = record.get("coordinate_validation", "pipeline crystallographic reference").replace("_", " ")
        self.reference_label.setText(f"Reference: {record.get('identity', 'not attached')} ({origin}).\n"
            "Saved reports retain their original reference; analyze again after changing it.")

    def attach_reference(self):
        if self.worker or self.reader or self.reference_worker:
            return
        filename, _ = QFileDialog.getOpenFileName(self, "Load co-crystallized ligand in the receptor coordinate frame", filter="SDF files (*.sdf)")
        if not filename:
            return
        answer = QMessageBox.question(self, "Crystallographic reference",
            "Confirm that this SDF contains one co-crystallized ligand pose, already in the prepared receptor's coordinate frame.\n\n"
            "Do not use a redocked or independently generated conformer. Coordinates will not be aligned or optimized. "
            "If the original crystal source is unavailable, this declaration will be recorded as user-supplied coordinate evidence.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel, QMessageBox.StandardButton.Cancel)
        if answer != QMessageBox.StandardButton.Yes:
            return
        from ..services.interaction_analysis import attach_crystallographic_reference
        self.analyze.setEnabled(False); self.load_reference.setEnabled(False)
        self.reference_worker = TaskWorker(lambda: attach_crystallographic_reference(
            self.repository, self.profile_id, filename, confirm_same_frame=True), self)
        self.reference_worker.succeeded.connect(self.reference_attached)
        self.reference_worker.failed.connect(self.failure)
        self.reference_worker.finished.connect(self.reference_finished)
        self.reference_worker.start()

    def reference_attached(self, record):
        self._last_error = None
        self.refresh_reference_label()
        self.status.setText("Co-crystallized reference attached. Run Analyze Interactions to generate a new comparison.")

    def reference_finished(self):
        worker, self.reference_worker = self.reference_worker, None
        worker.deleteLater()
        self.analyze.setEnabled(True); self.load_reference.setEnabled(True)
        if self._closing:
            self.finish_close()

    def start(self):
        if self.worker or self.reader or self.reference_worker:
            return
        self.analyze.setEnabled(False); self.scope.setEnabled(False); self.cancel.setEnabled(True)
        self.load_reference.setEnabled(False)
        self._last_error = None
        self.status.setText("Validating existing receptor, reference, and SDF exports…")
        self.worker = InteractionWorker(self.repository, self.profile_id, self.scope.currentData(), self.force.isChecked(), self)
        self.worker.progress.connect(lambda event: self.status.setText(event.message or f"Interaction analysis: {event.event}"))
        self.worker.succeeded.connect(lambda result: self.status.setText(f"{result[0]} poses analyzed; {result[1]} failed"))
        self.worker.failed.connect(self.failure)
        self.worker.finished.connect(self.finished_analysis)
        self.worker.start()

    def failure(self, message):
        self._last_error = message
        self.status.setText(message)

    def finished_analysis(self):
        worker, self.worker = self.worker, None
        worker.deleteLater()
        self.cancel.setEnabled(False); self.scope.setEnabled(True)
        self.load_reference.setEnabled(True)
        if self._closing:
            self.finish_close()
        else:
            self.load()

    def sort(self, column):
        key = self.COLUMNS[column][0]
        self._descending = not self._descending if key == self._sort else False
        self._sort, self._offset = key, 0
        self.table.horizontalHeader().setSortIndicator(column,
            Qt.SortOrder.DescendingOrder if self._descending else Qt.SortOrder.AscendingOrder)
        self.load()

    def page(self, amount):
        self._offset = max(0, self._offset + amount)
        self.load()

    def load(self):
        if self.reader:
            self._read_pending = True
            return
        if self.worker or self.reference_worker or self._closing:
            return
        self.analyze.setEnabled(False)
        self.previous.setEnabled(False); self.next.setEnabled(False)
        for button in self.downloads:
            button.setEnabled(False)
        scope, offset, sort, descending = self.scope.currentData(), self._offset, self._sort, self._descending
        def read():
            with self.repository.connection() as conn:
                snapshot = conn.execute("""SELECT a.* FROM interaction_runs a JOIN workflow_runs w
                    ON w.workflow_run_id=a.workflow_run_id WHERE a.receptor_profile_id=?
                    AND a.status IN ('completed','completed_with_errors')
                    AND (json_extract(w.settings_json,'$.requested_run_id') IS ?
                         OR (? IS NOT NULL AND json_extract(w.settings_json,'$.requested_run_id') IS NULL))
                    ORDER BY a.finished_at DESC,a.rowid DESC LIMIT 1""", (self.profile_id, scope, scope)).fetchone()
                if snapshot is None:
                    return None, [], 0
                condition = "analysis_id=?"
                parameters = [snapshot["analysis_id"]]
                if scope:
                    condition += " AND json_extract(summary_json,'$.run_id')=?"
                    parameters.append(scope)
                count = conn.execute(f"SELECT COUNT(*) FROM interaction_pose_summaries WHERE {condition}", parameters).fetchone()[0]
                # Sort keys come exclusively from COLUMNS, never free-form input.
                ordering = "DESC" if descending else "ASC"
                rows = conn.execute(f"SELECT summary_json FROM interaction_pose_summaries WHERE {condition} ORDER BY json_extract(summary_json,'$.{sort}') IS NULL, json_extract(summary_json,'$.{sort}') {ordering}, source_record_id LIMIT 200 OFFSET ?",
                                    parameters + [offset]).fetchall()
            return dict(snapshot), [json.loads(row[0]) for row in rows], count
        self.reader = TaskWorker(read, self)
        self.reader.succeeded.connect(self.loaded)
        self.reader.failed.connect(self.failure)
        self.reader.finished.connect(self.finished_read)
        self.reader.start()

    def loaded(self, result):
        snapshot, rows, count = result
        self.table.setRowCount(len(rows))
        for i, row in enumerate(rows):
            for j, (key, _) in enumerate(self.COLUMNS):
                value = row.get(key)
                item = QTableWidgetItem()
                item.setData(Qt.ItemDataRole.DisplayRole, "NA" if value is None else value)
                self.table.setItem(i, j, item)
        self._output = self.repository.root / snapshot["output_path"] if snapshot else None
        if snapshot:
            self.status.setText(f"Snapshot {snapshot['finished_at']}: {snapshot['succeeded']} poses analyzed; {snapshot['failed']} failed. "
                                f"Showing {self._offset + 1 if count else 0}–{min(self._offset + 200, count)} of {count}.")
            if self._last_error:
                self.status.setText(f"Latest attempt: {self._last_error}. Previous " + self.status.text())
        elif not self.status.text() or self.status.text().startswith("Loading"):
            self.status.setText("No saved analysis. Export SDF poses, then Analyze Interactions.")
        self.previous.setEnabled(self._offset > 0)
        self.next.setEnabled(self._offset + 200 < count)
        for button in self.downloads:
            button.setEnabled(self._output is not None)

    def finished_read(self):
        worker, self.reader = self.reader, None
        worker.deleteLater()
        self.analyze.setEnabled(True)
        if self._closing:
            self.finish_close()
        elif self._read_pending:
            self._read_pending = False
            self.load()

    def download(self, name):
        if self._output is None:
            return
        destination, _ = QFileDialog.getSaveFileName(self, "Save interaction report", name, "CSV files (*.csv)")
        if destination:
            try:
                shutil.copyfile(self._output / name, destination)
            except OSError as exc:
                QMessageBox.warning(self, "Report download failed", str(exc))

    def done(self, result):
        self._closing, self._close_result = True, result
        if self.worker:
            self.worker.requestInterruption()
        self.finish_close()

    def finish_close(self):
        if self.worker is None and self.reader is None and self.reference_worker is None:
            super().done(self._close_result)

    def closeEvent(self, event):
        self.done(QDialog.DialogCode.Rejected)
        if self.worker or self.reader or self.reference_worker:
            event.ignore()
        else:
            event.accept()
